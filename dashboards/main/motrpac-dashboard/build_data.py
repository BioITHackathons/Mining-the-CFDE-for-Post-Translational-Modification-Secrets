"""
Build motrpac-dashboard/data.js from MoTrPAC MotrpacRatTraining6moData R package.

Includes PTM sites (PHOSPHO, ACETYL, UBIQ) AND protein abundance (PROT).
Grouped by (feature_ID, tissue) for compact inline embedding.

Usage:
    uv run --with pyreadr --with requests --with pandas python build_data.py
"""

import json, os, re, tempfile, time
from collections import defaultdict
import pandas as pd
import pyreadr, requests

DATA_BASE = "https://github.com/MoTrPAC/MotrpacRatTraining6moData/raw/main/data"
OUT = os.path.join(os.path.dirname(__file__), "data.js")
SITE_RE = re.compile(r"^(.+?)_([A-Z])(\d+)")


def download_rda(filename):
    url = f"{DATA_BASE}/{filename}"
    print(f"  Downloading {filename}...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    tmp = tempfile.NamedTemporaryFile(suffix=".rda", delete=False)
    tmp.write(r.content)
    tmp.close()
    try:
        return pyreadr.read_r(tmp.name)
    finally:
        os.unlink(tmp.name)


def resolve_uniprot_batch(gene_symbols):
    gene_to_acc = {}
    batch_size = 80
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i + batch_size]
        query = " OR ".join(f"gene_exact:{g}" for g in batch)
        url = (
            f"https://rest.uniprot.org/uniprotkb/search"
            f"?query=({query})+AND+organism_id:9606+AND+reviewed:true"
            f"&fields=accession,gene_names&format=json&size=500"
        )
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                for entry in resp.json().get("results", []):
                    acc = entry.get("primaryAccession", "")
                    for gn in entry.get("genes", []):
                        name = gn.get("geneName", {}).get("value", "").upper()
                        if name and name not in gene_to_acc:
                            gene_to_acc[name] = acc
            time.sleep(0.3)
        except Exception as e:
            print(f"  Warning: batch failed: {e}")
    return gene_to_acc


def main():
    print("Downloading MoTrPAC data...")
    tr = download_rda("TRAINING_REGULATED_FEATURES.rda")["TRAINING_REGULATED_FEATURES"]
    f2g = download_rda("FEATURE_TO_GENE_FILT.rda")["FEATURE_TO_GENE_FILT"]
    r2h_gene = download_rda("RAT_TO_HUMAN_GENE.rda")["RAT_TO_HUMAN_GENE"]
    r2h_phospho = download_rda("RAT_TO_HUMAN_PHOSPHO.rda")["RAT_TO_HUMAN_PHOSPHO"]

    # ── Lookup maps ──────────────────────────────────────────────────────
    rat2human = {}
    for _, row in r2h_gene.iterrows():
        rs, hs = row.get("RAT_SYMBOL"), row.get("HUMAN_ORTHOLOG_SYMBOL")
        if rs and hs and str(hs) != "nan":
            rat2human[str(rs).upper()] = str(hs).upper()

    fid2gene = {}
    for _, row in f2g.iterrows():
        fid, gene = row.get("feature_ID"), row.get("gene_symbol")
        if fid and gene:
            fid2gene[str(fid)] = str(gene)

    rat2human_site = {}
    for _, row in r2h_phospho.iterrows():
        rid, hid = str(row.get("ptm_id_rat_refseq", "")), str(row.get("ptm_id_human_uniprot", ""))
        if rid and hid and hid != "nan":
            rat2human_site[rid] = hid

    # ── Resolve UniProt ──────────────────────────────────────────────────
    print("Resolving UniProt accessions...")
    all_human_genes = set()
    ptm_and_prot = tr[tr["assay"].isin(["PHOSPHO", "ACETYL", "UBIQ", "PROT"])]
    for _, row in ptm_and_prot.iterrows():
        rg = fid2gene.get(str(row["feature_ID"]), "")
        if rg:
            hg = rat2human.get(rg.upper(), "")
            if hg and hg != "nan":
                all_human_genes.add(hg)

    gene_to_uniprot = resolve_uniprot_batch(sorted(all_human_genes - {"nan", ""}))
    print(f"  Resolved: {len(gene_to_uniprot)}")

    for fid, hid in rat2human_site.items():
        parts = hid.split("_", 1)
        if len(parts) == 2:
            rg = fid2gene.get(fid, "")
            if rg:
                hg = rat2human.get(rg.upper(), "")
                if hg and hg not in gene_to_uniprot:
                    gene_to_uniprot[hg] = parts[0]

    # ── Process PTM rows ─────────────────────────────────────────────────
    ptm = tr[tr["assay"].isin(["PHOSPHO", "ACETYL", "UBIQ"])].copy()
    print(f"PTM rows: {len(ptm)}")

    records = []
    for _, row in ptm.iterrows():
        fid = str(row["feature_ID"])
        assay, tissue, sex, time_pt = row["assay"], row["tissue"], row["sex"], row["training_group"]
        logfc, pval = row["timewise_logFC"], row["timewise_p_value"]

        rg = fid2gene.get(fid, "")
        hg = rat2human.get(rg.upper(), rg) if rg else ""

        # Always parse the rat site from the feature_ID
        rat_residue, rat_position = "", 0
        m_rat = SITE_RE.match(fid)
        if m_rat:
            rat_residue, rat_position = m_rat.group(2), int(m_rat.group(3))

        if rat_position == 0:
            continue

        # Try to get confirmed human mapping
        human_residue, human_position, mapped = "", 0, 0
        uniprot = ""

        if assay == "PHOSPHO":
            hid = rat2human_site.get(fid, "")
            if hid:
                parts = hid.split("_", 1)
                if len(parts) == 2:
                    uniprot = parts[0]
                    m = re.match(r"([A-Z])(\d+)", parts[1])
                    if m:
                        human_residue = m.group(1)
                        human_position = int(m.group(2))
                        mapped = 1

        # If no confirmed human mapping, use rat coordinates as estimate
        if human_position == 0:
            human_residue = rat_residue
            human_position = rat_position

        if not uniprot and hg:
            uniprot = gene_to_uniprot.get(hg.upper(), "")

        records.append(dict(fid=fid, assay=assay, tissue=tissue, sex=sex, time=time_pt,
                            logFC=round(float(logfc), 3) if pd.notna(logfc) else None,
                            pval=round(float(pval), 4) if pd.notna(pval) else None,
                            gene=hg, residue=human_residue, position=human_position,
                            ratRes=rat_residue, ratPos=rat_position,
                            mapped=mapped, uniprot=uniprot))

    # ── Process PROT (protein abundance) rows ────────────────────────────
    prot = tr[tr["assay"] == "PROT"].copy()
    print(f"PROT rows: {len(prot)}")

    prot_records = []
    for _, row in prot.iterrows():
        fid = str(row["feature_ID"])
        tissue, sex, time_pt = row["tissue"], row["sex"], row["training_group"]
        logfc, pval = row["timewise_logFC"], row["timewise_p_value"]

        rg = fid2gene.get(fid, "")
        hg = rat2human.get(rg.upper(), rg) if rg else ""
        uniprot = gene_to_uniprot.get(hg.upper(), "") if hg else ""

        if not hg or hg == "nan":
            continue

        prot_records.append(dict(gene=hg, tissue=tissue, sex=sex, time=time_pt,
                                 uniprot=uniprot,
                                 logFC=round(float(logfc), 3) if pd.notna(logfc) else None,
                                 pval=round(float(pval), 4) if pd.notna(pval) else None))

    print(f"PROT records (mapped): {len(prot_records)}")

    # ── Group PTM data ───────────────────────────────────────────────────
    groups = defaultdict(dict)
    group_meta = {}
    for r in records:
        key = (r["fid"], r["tissue"])
        ck = f"{r['sex'][0]}{r['time']}"
        groups[key][ck] = [r["logFC"], r["pval"]]
        if key not in group_meta:
            group_meta[key] = r

    grouped = []
    mapped_count = 0
    for key, conds in groups.items():
        ref = group_meta[key]
        entry = {"f": ref["fid"], "t": ref["tissue"], "a": ref["assay"][0],
                 "g": ref["gene"], "r": ref["residue"], "P": ref["position"],
                 "u": ref["uniprot"], "v": conds,
                 "rR": ref["ratRes"], "rP": ref["ratPos"], "m": ref["mapped"]}
        grouped.append(entry)
        if ref["mapped"]:
            mapped_count += 1

    print(f"Grouped PTM sites: {len(grouped)} ({mapped_count} with confirmed human mapping)")

    # ── Group PROT data by (gene, tissue) ────────────────────────────────
    prot_groups = defaultdict(dict)
    prot_meta = {}
    for r in prot_records:
        key = (r["gene"], r["tissue"])
        ck = f"{r['sex'][0]}{r['time']}"
        prot_groups[key][ck] = [r["logFC"], r["pval"]]
        if key not in prot_meta:
            prot_meta[key] = r

    prot_grouped = []
    for key, conds in prot_groups.items():
        ref = prot_meta[key]
        prot_grouped.append({"g": ref["gene"], "t": ref["tissue"],
                              "u": ref["uniprot"], "v": conds})

    print(f"Grouped PROT entries: {len(prot_grouped)}")

    # ── Protein index ────────────────────────────────────────────────────
    proteins = {}
    for g in grouped:
        pkey = g["u"] or g["g"]
        if not pkey or pkey == "nan":
            continue
        if pkey not in proteins:
            proteins[pkey] = {"uniprot": g["u"], "gene": g["g"],
                              "sites": set(), "tissues": set(), "assays": set()}
        proteins[pkey]["sites"].add(g["P"])
        proteins[pkey]["tissues"].add(g["t"])
        proteins[pkey]["assays"].add({"P": "PHOSPHO", "A": "ACETYL", "U": "UBIQ"}[g["a"]])

    # Mark which proteins also have abundance data
    prot_gene_set = set(r["gene"] for r in prot_records if r["gene"])
    protein_index = []
    for pkey, p in sorted(proteins.items(), key=lambda x: len(x[1]["sites"]), reverse=True):
        protein_index.append({
            "id": pkey, "uniprot": p["uniprot"], "gene": p["gene"],
            "nSites": len(p["sites"]), "tissues": sorted(p["tissues"]),
            "assays": sorted(p["assays"]),
            "hasProt": p["gene"] in prot_gene_set,
        })

    has_prot_count = sum(1 for p in protein_index if p["hasProt"])
    print(f"Proteins: {len(protein_index)} ({has_prot_count} with abundance data)")

    # ── Metadata ─────────────────────────────────────────────────────────
    meta = {
        "tissues": sorted(set(g["t"] for g in grouped)),
        "assays": ["ACETYL", "PHOSPHO", "UBIQ"],
        "times": ["1w", "2w", "4w", "8w"],
        "sexes": ["female", "male"],
        "protein_count": len(protein_index),
        "record_count": len(records),
        "tissue_labels": {
            "CORTEX": "Cortex", "HEART": "Heart", "KIDNEY": "Kidney",
            "LIVER": "Liver", "LUNG": "Lung", "SKM-GN": "Gastrocnemius",
            "WAT-SC": "Subcutaneous WAT",
        },
        "assay_labels": {
            "PHOSPHO": "Phosphorylation", "ACETYL": "Acetylation",
            "UBIQ": "Ubiquitylation",
        },
    }

    # ── Write JS ─────────────────────────────────────────────────────────
    sep = (",", ":")
    js = "// Auto-generated by build_data.py\n"
    js += f"const MOTRPAC_META = {json.dumps(meta, indent=2)};\n\n"
    js += f"const MOTRPAC_PROTEINS = {json.dumps(protein_index, separators=sep)};\n\n"
    js += f"const MOTRPAC_GROUPED = {json.dumps(grouped, separators=sep)};\n\n"
    js += f"const MOTRPAC_PROT = {json.dumps(prot_grouped, separators=sep)};\n"

    with open(OUT, "w") as f:
        f.write(js)
    print(f"Wrote {OUT} ({len(js):,} bytes)")


if __name__ == "__main__":
    main()
