#!/usr/bin/env python3
"""
Resumable, chunk-based database builder.

Runs each stage independently. If a stage fails partway through, re-running
the script skips already-processed items. Progress is tracked in the DB itself
and in checkpoint files under .checkpoints/.

Usage:
    python build_data.py                    # run all stages
    python build_data.py --stage seed       # just seed
    python build_data.py --stage glygen     # just GlyGen (resumable)
    python build_data.py --stage motrpac    # just MoTrPAC
    python build_data.py --stage kc         # just CFDE KC
    python build_data.py --stage summary    # just print summary
"""

import argparse
import gc
import json
import os
import re
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

import requests

# Keep memory footprint small: disable requests connection pooling for long runs
SESSION = requests.Session()
SESSION.headers.update({"Connection": "close"})

DB_PATH = Path(__file__).resolve().parent / "ptm_disease.db"
SCHEMA_PATH = Path(__file__).resolve().parent / "ptm_disease_db_schema.sql"
CHECKPOINT_DIR = Path(__file__).resolve().parent / ".checkpoints"

sys.path.insert(0, str(Path(__file__).resolve().parent / "loaders"))
from db_utils import *
from seed_proteins import SEED_PROTEINS, SEED_DISEASES


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def get_db():
    """Get or create the database."""
    if not DB_PATH.exists():
        print(f"Creating database at {DB_PATH}...")
        conn = sqlite3.connect(str(DB_PATH))
        conn.executescript(SCHEMA_PATH.read_text())
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA foreign_keys = ON")
        conn.commit()
        conn.close()
    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA foreign_keys = ON")
    conn.row_factory = sqlite3.Row
    return conn


def load_checkpoint(name):
    """Load a set of completed items from a checkpoint file."""
    path = CHECKPOINT_DIR / f"{name}.json"
    if path.exists():
        return set(json.loads(path.read_text()))
    return set()


def save_checkpoint(name, done_set):
    """Save completed items to a checkpoint file."""
    CHECKPOINT_DIR.mkdir(exist_ok=True)
    path = CHECKPOINT_DIR / f"{name}.json"
    path.write_text(json.dumps(sorted(done_set)))


AA = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
      "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
      "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}

GLYCO_MAP = {
    "N-linked": ("MOD:00006", "N-linked glycosylation", "glycosylation"),
    "O-linked": ("MOD:00007", "O-linked glycosylation", "glycosylation"),
    "C-linked": ("MOD:01084", "C-linked glycosylation", "glycosylation"),
}


def extract_pmids(ev_list):
    return [e["id"] for e in ev_list if e.get("database") == "PubMed"]


def ev_type(ev_list):
    dbs = {e.get("database", "").lower() for e in ev_list}
    if "pubmed" in dbs: return "experimental"
    if "uniprotkb" in dbs or "iptmnet" in dbs: return "curated"
    return "inferred"


# ═══════════════════════════════════════════════════════════════════════════
# Stage 1: Seed
# ═══════════════════════════════════════════════════════════════════════════

def stage_seed():
    print("\n" + "=" * 60)
    print("STAGE: Seed reference data")
    print("=" * 60)
    conn = get_db()

    n_existing = conn.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
    if n_existing >= len(SEED_PROTEINS):
        print(f"  Already seeded ({n_existing} proteins). Skipping.")
        conn.close()
        return

    for acc, gene in SEED_PROTEINS.items():
        upsert_protein(conn, acc, gene)
    for doid, info in SEED_DISEASES.items():
        upsert_disease(conn, doid, info["name"], category=info["category"])
        for acc in info["proteins"]:
            upsert_protein_disease(conn, acc, doid, association_type="genetic")
    conn.commit()
    print(f"  Seeded {len(SEED_PROTEINS)} proteins, {len(SEED_DISEASES)} diseases.")
    conn.close()


# ═══════════════════════════════════════════════════════════════════════════
# Stage 2: GlyGen (resumable, chunked)
# ═══════════════════════════════════════════════════════════════════════════

def glygen_discover_accessions():
    """Discover all human protein accessions with glycosylation data."""
    cache = CHECKPOINT_DIR / "glygen_accessions.json"
    if cache.exists():
        accs = json.loads(cache.read_text())
        print(f"  Loaded {len(accs)} cached accessions from {cache}")
        return accs

    print("  Discovering human glycoproteins from GlyGen search API...")
    r = requests.post("https://api.glygen.org/protein/search",
                      json={"organism": {"id": 9606}, "glycosylation_evidence": "reported"}, timeout=30)
    r.raise_for_status()
    list_id = r.json()["list_id"]

    r2 = requests.post("https://api.glygen.org/protein/list",
                        json={"id": list_id, "offset": 1, "limit": 1}, timeout=30)
    total = r2.json()["pagination"]["total_length"]
    print(f"  Total in GlyGen: {total}")

    accessions = set()
    offset = 1
    while offset <= total:
        r3 = requests.post("https://api.glygen.org/protein/list",
                           json={"id": list_id, "offset": offset, "limit": 500}, timeout=60)
        for p in r3.json().get("results", []):
            acc = p.get("uniprot_canonical_ac", "").split("-")[0]
            n = (p.get("reported_n_glycosites") or 0) + (p.get("reported_o_glycosites") or 0)
            if n > 0 and acc:
                accessions.add(acc)
        offset += 500
        time.sleep(0.2)

    result = sorted(accessions)
    CHECKPOINT_DIR.mkdir(exist_ok=True)
    cache.write_text(json.dumps(result))
    print(f"  Found {len(result)} proteins with glycosylation sites. Cached.")
    return result


def glygen_load_one(conn, acc):
    """Load all PTM data for one protein from GlyGen detail API."""
    try:
        r = requests.post(f"https://api.glygen.org/protein/detail/{acc}", json={}, timeout=60)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"    WARN {acc}: {e}")
        return False

    gene_names = data.get("gene_names", [])
    gene = gene_names[0].get("name", acc) if gene_names else acc
    pnames = data.get("protein_names", [])
    pname = pnames[0].get("name") if pnames else None
    sl = len(data.get("sequence", {}).get("sequence", "")) or None
    upsert_protein(conn, acc, gene, pname, sequence_length=sl)

    # Diseases
    for d in data.get("disease", []):
        doid = d.get("disease_id", "")
        rn = d.get("recommended_name", {})
        if doid.startswith("DOID:"):
            upsert_disease(conn, doid, rn.get("name", doid), rn.get("description"))
            upsert_protein_disease(conn, acc, doid, source_id="glygen")

    # Glycosylation
    upsert_ptm_type(conn, "MOD:00696", "phosphorylation", "phosphorylation")
    for site in data.get("glycosylation", []):
        gt = site.get("type", "O-linked")
        pid, pnm, pc = GLYCO_MAP.get(gt, ("MOD:00007", gt, "glycosylation"))
        upsert_ptm_type(conn, pid, pnm, pc)
        pos = site.get("start_pos")
        res = AA.get(site.get("residue") or site.get("start_aa"), site.get("residue"))
        if not pos:
            continue
        sid = upsert_ptm_site(conn, acc, pos, res, pid,
                              subtype=site.get("subtype"), flanking_sequence=site.get("site_seq"))
        evl = site.get("evidence", [])
        pms = extract_pmids(evl)
        et = ev_type(evl)
        # Store only first PubMed ID as evidence (skip raw_record to save memory/disk)
        if pms:
            insert_evidence(conn, sid, "glygen", et, pubmed_id=pms[0])
        else:
            insert_evidence(conn, sid, "glygen", et)
        gtc = site.get("glytoucan_ac")
        if gtc:
            upsert_modifier(conn, gtc, "glycan")
            upsert_ptm_site_modifier(conn, sid, gtc, source_id="glygen", evidence_type=et,
                                     pubmed_id=pms[0] if pms else None)

    # Phosphorylation
    for site in data.get("phosphorylation", []):
        pos = site.get("start_pos")
        res = AA.get(site.get("residue"), site.get("residue"))
        if not pos:
            continue
        sid = upsert_ptm_site(conn, acc, pos, res, "MOD:00696", notes=site.get("comment"))
        evl = site.get("evidence", [])
        pms = extract_pmids(evl)
        et = ev_type(evl)
        if pms:
            insert_evidence(conn, sid, "glygen", et, pubmed_id=pms[0])
        else:
            insert_evidence(conn, sid, "glygen", et)

    # SNVs
    for snv in data.get("snv", []):
        pos = snv.get("start_pos")
        if not pos:
            continue
        dbsnp = None
        for ev in snv.get("evidence", []):
            if ev.get("database") == "dbSNP":
                dbsnp = ev.get("id"); break
        vid = upsert_variant(conn, acc, pos, ref_residue=snv.get("sequence_org"),
                             alt_residue=snv.get("sequence_mut"), dbsnp_id=dbsnp, source_id="glygen")
        for row in conn.execute("SELECT ptm_site_id FROM ptm_site WHERE uniprot_accession=? AND position=?",
                                (acc, pos)).fetchall():
            eff = (snv.get("glycoeffect", []) or ["unknown"])[0]
            conn.execute("INSERT OR IGNORE INTO variant_ptm_site (variant_id,ptm_site_id,effect,distance,notes) VALUES(?,?,?,0,?)",
                         (vid, row[0], eff, snv.get("comment")))
        for d in snv.get("disease", []):
            doid = d.get("disease_id", "")
            if doid.startswith("DOID:"):
                rn = d.get("recommended_name", {})
                upsert_disease(conn, doid, rn.get("name", doid), rn.get("description"))
    return True


def stage_glygen():
    print("\n" + "=" * 60)
    print("STAGE: GlyGen (all human glycoproteins)")
    print("=" * 60)

    accessions = glygen_discover_accessions()
    done = load_checkpoint("glygen_done")
    remaining = [a for a in accessions if a not in done]

    print(f"  Total: {len(accessions)}, Done: {len(done)}, Remaining: {len(remaining)}")
    if not remaining:
        print("  All done. Skipping.")
        return

    conn = get_db()
    upsert_data_source(conn, "glygen", "GlyGen", url="https://glygen.org",
                       cfde_program="GlyGen", description="Glycosylation and PTM data from GlyGen")
    conn.commit()

    BATCH = 50
    for i in range(0, len(remaining), BATCH):
        batch = remaining[i:i+BATCH]
        batch_num = i // BATCH + 1
        total_batches = (len(remaining) + BATCH - 1) // BATCH
        print(f"\n  Batch {batch_num}/{total_batches} ({i+1}-{i+len(batch)} of {len(remaining)})...")
        for acc in batch:
            glygen_load_one(conn, acc)
            done.add(acc)
            time.sleep(0.15)  # lighter rate limit
        conn.commit()
        save_checkpoint("glygen_done", done)
        n_sites = conn.execute("SELECT COUNT(*) FROM ptm_site").fetchone()[0]
        n_proteins = conn.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
        print(f"    Committed. DB: {n_proteins} proteins, {n_sites} PTM sites")
        gc.collect()  # free any dangling JSON responses

    conn.close()
    print(f"\n  GlyGen complete. Processed {len(accessions)} proteins.")


# ═══════════════════════════════════════════════════════════════════════════
# Stage 3: MoTrPAC
# ═══════════════════════════════════════════════════════════════════════════

def stage_motrpac():
    print("\n" + "=" * 60)
    print("STAGE: MoTrPAC (phospho, acetyl, ubiquityl)")
    print("=" * 60)

    done = load_checkpoint("motrpac_done")
    if done:
        print(f"  Already completed ({len(done)} sites). Skipping.")
        return

    try:
        import pyreadr
    except ImportError:
        print("  ERROR: pyreadr not installed. Run: uv pip install pyreadr")
        return

    conn = get_db()
    upsert_data_source(conn, "motrpac", "MoTrPAC", url="https://motrpac-data.org",
                       cfde_program="MoTrPAC",
                       description="PTM data from MoTrPAC endurance training study. Nature 629, 174-183 (2024).")
    for pid, pnm, pc in [("MOD:00696","phosphorylation","phosphorylation"),
                          ("MOD:00394","acetylation","acetylation"),
                          ("MOD:01148","ubiquitylation","ubiquitylation")]:
        upsert_ptm_type(conn, pid, pnm, pc)
    conn.commit()

    BASE = "https://github.com/MoTrPAC/MotrpacRatTraining6moData/raw/main/data"
    ASSAY_MAP = {"PHOSPHO":("MOD:00696","phosphorylation"),
                 "ACETYL":("MOD:00394","acetylation"),
                 "UBIQ":("MOD:01148","ubiquitylation")}
    SITE_RE = re.compile(r"^(.+?)_([A-Z])(\d+)([a-z])(.*)$")

    def dl(fn):
        cache = CHECKPOINT_DIR / fn
        if cache.exists():
            print(f"  Using cached {fn}")
            return pyreadr.read_r(str(cache))
        print(f"  Downloading {fn}...")
        r = requests.get(f"{BASE}/{fn}", timeout=120); r.raise_for_status()
        CHECKPOINT_DIR.mkdir(exist_ok=True)
        cache.write_bytes(r.content)
        return pyreadr.read_r(str(cache))

    # Load and immediately extract only what we need, then free the DataFrames
    print("  Loading rat-to-human gene map...")
    r2h_g = dl("RAT_TO_HUMAN_GENE.rda")["RAT_TO_HUMAN_GENE"]
    rat2hum = dict(zip(r2h_g["RAT_SYMBOL"].str.upper(), r2h_g["HUMAN_ORTHOLOG_SYMBOL"].str.upper()))
    rat2hum = {k: v for k, v in rat2hum.items() if v and v != "NAN"}
    del r2h_g; gc.collect()

    print("  Loading rat-to-human phospho site map...")
    r2h_p = dl("RAT_TO_HUMAN_PHOSPHO.rda")["RAT_TO_HUMAN_PHOSPHO"]
    rat2hum_ptm = dict(zip(r2h_p["ptm_id_rat_refseq"].astype(str),
                           r2h_p["ptm_id_human_uniprot"].astype(str)))
    rat2hum_ptm = {k: v for k, v in rat2hum_ptm.items() if v and v != "nan"}
    del r2h_p; gc.collect()

    print("  Loading feature-to-gene map...")
    f2g = dl("FEATURE_TO_GENE_FILT.rda")["FEATURE_TO_GENE_FILT"]
    feat2gene = dict(zip(f2g["feature_ID"].astype(str), f2g["gene_symbol"].astype(str)))
    del f2g; gc.collect()

    print("  Loading training-regulated features (PTM only)...")
    tr_full = dl("TRAINING_REGULATED_FEATURES.rda")["TRAINING_REGULATED_FEATURES"]
    # Filter to PTM assays immediately and keep only needed columns
    ptm = tr_full[tr_full["assay"].isin(["PHOSPHO","ACETYL","UBIQ"])][
        ["feature_ID","assay","tissue","timewise_p_value"]
    ].copy()
    del tr_full; gc.collect()

    # Pre-compute min p-value per (feature, tissue)
    min_p = ptm.groupby(["feature_ID","tissue"])["timewise_p_value"].min().to_dict()
    uniq = ptm.drop_duplicates(subset=["feature_ID","assay","tissue"])
    print(f"  Unique PTM feature/tissue combos: {len(uniq)}")
    del ptm; gc.collect()

    # Build gene->acc map from current DB
    gene2acc = {g.upper(): a for a, g in SEED_PROTEINS.items()}
    for row in conn.execute("SELECT uniprot_accession, gene_symbol FROM protein"):
        gene2acc.setdefault(row[1].upper() if row[1] else "", row[0])

    loaded = 0
    seen = set()
    for fid, assay, tissue in zip(uniq["feature_ID"].astype(str),
                                   uniq["assay"], uniq["tissue"]):
        rg = feat2gene.get(fid)
        if not rg or not isinstance(rg, str): continue
        hg = rat2hum.get(rg.upper())
        if not hg or not isinstance(hg, str): continue
        acc = gene2acc.get(hg.upper())
        if not acc: continue
        pid, pnm = ASSAY_MAP[assay]
        pos = res = None
        if assay == "PHOSPHO":
            hp = rat2hum_ptm.get(fid)
            if hp and isinstance(hp, str):
                parts = hp.split("_", 1)
                if len(parts) == 2:
                    m = re.match(r"([A-Z])(\d+)", parts[1])
                    if m: res, pos = m.group(1), int(m.group(2))
        if pos is None:
            m = SITE_RE.match(fid)
            if m: res, pos = m.group(2), int(m.group(3))
            else: continue
        sk = (acc, pos, pid)
        if sk in seen: continue
        seen.add(sk)
        upsert_protein(conn, acc, hg)
        sid = upsert_ptm_site(conn, acc, pos, res, pid, notes=f"MoTrPAC; rat={rg}")
        mp = min_p.get((fid, tissue))
        insert_evidence(conn, sid, "motrpac", "experimental", method="mass spectrometry (TMT)",
                        confidence_score=1-mp if mp is not None and mp==mp else None,
                        raw_record={"feature_ID":fid,"assay":assay,"tissue":tissue,
                                    "rat_gene":rg,"human_gene":hg})
        loaded += 1

    del uniq, min_p, feat2gene, rat2hum, rat2hum_ptm; gc.collect()
    conn.commit()
    save_checkpoint("motrpac_done", list(seen))
    print(f"  MoTrPAC complete: {loaded} PTM sites loaded.")
    conn.close()


# ═══════════════════════════════════════════════════════════════════════════
# Stage 4: CFDE KC (GTEx, Kids First, IDG) -- chunked
# ═══════════════════════════════════════════════════════════════════════════

def stage_kc():
    print("\n" + "=" * 60)
    print("STAGE: CFDE KC (GTEx, Kids First, IDG)")
    print("=" * 60)

    conn = get_db()
    for s in [("glygen_kc","GlyGen (KC)","GlyGen"),("gtex","GTEx","GTEx"),
              ("kids_first","Kids First","Kids First"),("idg","IDG / Pharos","IDG")]:
        upsert_data_source(conn, s[0], s[1], cfde_program=s[2])
    conn.commit()

    # All genes in DB
    genes = [r[0] for r in conn.execute("SELECT DISTINCT gene_symbol FROM protein").fetchall() if r[0]]
    done = load_checkpoint("kc_done")
    remaining = [g for g in genes if g not in done]
    print(f"  Genes: {len(genes)}, Done: {len(done)}, Remaining: {len(remaining)}")
    if not remaining:
        print("  All done. Skipping.")
        conn.close()
        return

    KC = "https://cfde.hugeampkpnbi.org/api/bio/query"
    GLYCO_KC = {"N-linked":("MOD:00006","N-linked glycosylation"),"O-linked":("MOD:00007","O-linked glycosylation")}

    def kc(idx, g):
        try:
            r = requests.get(f"{KC}/{idx}?q={g}", timeout=60); r.raise_for_status()
            return r.json().get("data", [])
        except: return []

    BATCH = 25
    for i in range(0, len(remaining), BATCH):
        batch = remaining[i:i+BATCH]
        print(f"\n  Batch {i//BATCH+1} ({i+1}-{i+len(batch)} of {len(remaining)})...")
        for gene in batch:
            # GTEx
            for row in kc("gtex-tstat", gene):
                upsert_tissue_expression(conn, gene, row.get("biosample","").replace("%27","'"),
                                         expression_tpm=row.get("tstat"), source_id="gtex")
            # Kids First
            ar = conn.execute("SELECT uniprot_accession FROM protein WHERE gene_symbol=? LIMIT 1", (gene,)).fetchone()
            if ar:
                acc = ar[0]
                for row in kc("kids-first-gene-variants", gene):
                    for vep in row.get("veprecords", []):
                        try: ppos = int(vep.get("Protein_position",""))
                        except: continue
                        aa = vep.get("Amino_acids",""); ra=al=None
                        if "/" in aa: ra,al = aa.split("/",1)
                        dbsnp = None
                        for hp in row.get("hprecords",[]):
                            if hp.get("id","").startswith("rs"): dbsnp=hp["id"]; break
                        vid = upsert_variant(conn, acc, ppos, ref_residue=ra, alt_residue=al,
                                             dbsnp_id=dbsnp or row.get("varId"),
                                             clinical_significance=vep.get("Consequence",""), source_id="kids_first")
                        for pr in conn.execute("SELECT ptm_site_id,position FROM ptm_site WHERE uniprot_accession=? AND ABS(position-?)<=5",
                                               (acc,ppos)).fetchall():
                            d=abs(ppos-pr[1])
                            conn.execute("INSERT OR IGNORE INTO variant_ptm_site VALUES(NULL,?,?,?,?,?)",
                                         (vid,pr[0],"abolishes_site" if d==0 else "adjacent",d,vep.get("Consequence","")))
                        break
            # IDG
            try:
                r = requests.post("https://pharos-api.ncats.io/graphql",
                                  json={"query":f'query{{target(q:{{sym:"{gene}"}}){{name tdl fam}}}}'}, timeout=30)
                r.raise_for_status()
                t = r.json().get("data",{}).get("target")
                if t and ar:
                    conn.execute("UPDATE protein SET protein_name=COALESCE(protein_name,?) WHERE uniprot_accession=?",
                                 (t.get("name"), ar[0]))
            except: pass

            done.add(gene)
            time.sleep(0.2)

        conn.commit()
        save_checkpoint("kc_done", done)
        print(f"    Committed. {len(done)}/{len(genes)} genes done.")

    conn.close()
    print(f"\n  KC complete.")


# ═══════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════

def stage_summary():
    print("\n" + "=" * 60)
    print("DATABASE SUMMARY")
    print("=" * 60)
    conn = get_db()
    for t, l in [("protein","Proteins"),("disease","Diseases"),("ptm_type","PTM types"),
                 ("data_source","Data sources"),("ptm_site","PTM sites"),("ptm_evidence","Evidence records"),
                 ("protein_disease","Protein-disease links"),("modifier","Modifiers"),
                 ("ptm_site_modifier","Site-modifier links"),("variant","Variants"),
                 ("variant_ptm_site","Variant-PTM overlaps"),("tissue_expression","Tissue expression")]:
        n = conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"  {l:.<40} {n:>8,}")

    print("\n  PTM sites by type:")
    for r in conn.execute("SELECT pt.name,COUNT(*) n FROM ptm_site ps JOIN ptm_type pt ON pt.ptm_type_id=ps.ptm_type_id GROUP BY pt.name ORDER BY n DESC"):
        print(f"    {r[0]:.<38} {r[1]:>8,}")

    print("\n  Top 15 proteins by PTM site count:")
    for r in conn.execute("SELECT p.gene_symbol, COUNT(*) n FROM ptm_site ps JOIN protein p ON p.uniprot_accession=ps.uniprot_accession GROUP BY p.gene_symbol ORDER BY n DESC LIMIT 15"):
        print(f"    {r[0]:.<38} {r[1]:>8,}")

    print(f"\n  Database size: {DB_PATH.stat().st_size / 1024 / 1024:.1f} MB")
    conn.close()


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

STAGES = {
    "seed": stage_seed,
    "glygen": stage_glygen,
    "motrpac": stage_motrpac,
    "kc": stage_kc,
    "summary": stage_summary,
}

def main():
    parser = argparse.ArgumentParser(description="Resumable PTM-disease DB builder")
    parser.add_argument("--stage", choices=list(STAGES.keys()), help="Run a single stage")
    parser.add_argument("--reset", action="store_true", help="Delete DB and checkpoints, start fresh")
    args = parser.parse_args()

    if args.reset:
        if DB_PATH.exists(): DB_PATH.unlink(); print(f"Deleted {DB_PATH}")
        if CHECKPOINT_DIR.exists():
            import shutil; shutil.rmtree(CHECKPOINT_DIR); print(f"Deleted {CHECKPOINT_DIR}")

    if args.stage:
        STAGES[args.stage]()
    else:
        for name, fn in STAGES.items():
            fn()

    print("\nDone.")


if __name__ == "__main__":
    main()
