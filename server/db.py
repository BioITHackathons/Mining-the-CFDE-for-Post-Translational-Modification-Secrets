"""
Live DB-backed query layer for the dashboard API.

Each function returns data in the same shape as the legacy JSON file it
replaces, so index.html doesn't need to change. NLP-derived fields (which
aren't in the DB yet) are merged in from glygen_data.json for the curated
proteins that have them; other proteins get empty lists for those fields.
"""

from __future__ import annotations

import json
import re
import sqlite3
from functools import lru_cache
from pathlib import Path

from server import glygen, pubmed

REPO_ROOT = Path(__file__).resolve().parents[1]
DB_PATH = REPO_ROOT / "pipeline" / "ptm_disease.db"
DASHBOARD_DIR = REPO_ROOT / "dashboards" / "main" / "ptm-dashboard"
NLP_AUG_PATH = DASHBOARD_DIR / "glygen_data.json"

NLP_FIELDS = (
    "expression_disease",
    "expression_tissue",
    "snv_with_disease_or_glycoeffect",
    "ptm_annotations_raw",
    "ptm_annotations_extracted",
    "publication_count",
)

# 1-letter -> 3-letter amino acid (matches the labels build_pqtl.py emitted).
AA_3 = {"A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","Q":"Gln","E":"Glu",
        "G":"Gly","H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe",
        "P":"Pro","S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val"}


def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA query_only = ON")
    return conn


@lru_cache(maxsize=1)
def _load_nlp_augmentation() -> dict:
    """Load NLP-derived fields keyed by uniprot_ac from the legacy JSON.

    Returns {} if the file isn't present.
    """
    if not NLP_AUG_PATH.exists():
        return {}
    raw = json.loads(NLP_AUG_PATH.read_text())
    aug = {}
    for p in raw.get("proteins", []):
        ac = p.get("uniprot_ac")
        if not ac:
            continue
        aug[ac] = {field: p.get(field, [] if field != "publication_count" else 0)
                   for field in NLP_FIELDS}
    return aug


def fetch_data() -> dict:
    """Return {summary, proteins[]} matching the shape of data.json.

    Protein universe comes from the live DB; NLP-derived fields are merged
    in from glygen_data.json where available, empty otherwise.
    """
    aug = _load_nlp_augmentation()
    conn = get_conn()

    # All proteins
    proteins_by_ac: dict[str, dict] = {}
    for row in conn.execute(
        "SELECT uniprot_accession, gene_symbol, protein_name, sequence_length "
        "FROM protein ORDER BY gene_symbol, uniprot_accession"
    ):
        ac = row["uniprot_accession"]
        gene = row["gene_symbol"] or ""
        pname = row["protein_name"] or ""
        display = f"{gene} ({pname})" if gene and pname else (gene or pname or ac)
        proteins_by_ac[ac] = {
            "uniprot_ac": ac,
            "display_name": display,
            "gene_name": gene,
            "protein_name": pname,
            "length": row["sequence_length"] or 0,
            "species": "Homo sapiens",
            "glycosylation_sites": [],
            "phosphorylation_sites": [],
            "diseases": [],
        }

    # Sites — single scan, route by ptm_type.category
    site_query = """
        SELECT ps.ptm_site_id, ps.uniprot_accession, ps.position, ps.residue,
               ps.subtype, ps.flanking_sequence, ps.notes,
               pt.name AS ptm_name, pt.category AS ptm_category
        FROM ptm_site ps
        JOIN ptm_type pt ON pt.ptm_type_id = ps.ptm_type_id
    """
    # First pass: collect site rows so we can attach glycan modifiers in one batch
    sites_rows = list(conn.execute(site_query))

    # Glycan modifier lookup: ptm_site_id -> first glytoucan_ac
    glycan_by_site: dict[int, str] = {}
    for row in conn.execute(
        "SELECT psm.ptm_site_id, m.modifier_id "
        "FROM ptm_site_modifier psm "
        "JOIN modifier m ON m.modifier_id = psm.modifier_id "
        "WHERE m.modifier_type = 'glycan'"
    ):
        glycan_by_site.setdefault(row["ptm_site_id"], row["modifier_id"])

    # Evidence count per site
    evcount_by_site: dict[int, int] = {}
    for row in conn.execute(
        "SELECT ptm_site_id, COUNT(*) AS n FROM ptm_evidence GROUP BY ptm_site_id"
    ):
        evcount_by_site[row["ptm_site_id"]] = row["n"]

    for s in sites_rows:
        ac = s["uniprot_accession"]
        if ac not in proteins_by_ac:
            continue
        sid = s["ptm_site_id"]
        pos = s["position"]
        res = s["residue"] or ""
        site_lbl = f"{res}{pos}" if res else str(pos)
        evc = evcount_by_site.get(sid, 0)
        if s["ptm_category"] == "phosphorylation":
            proteins_by_ac[ac]["phosphorylation_sites"].append({
                "start_pos": pos,
                "residue": res,
                "site_lbl": site_lbl,
                "comment": s["notes"] or "",
                "evidence_count": evc,
            })
        else:
            proteins_by_ac[ac]["glycosylation_sites"].append({
                "start_pos": pos,
                "end_pos": pos,
                "residue": res,
                "site_lbl": site_lbl,
                "type": s["ptm_name"] or "",
                "subtype": s["subtype"] or "",
                "glytoucan_ac": glycan_by_site.get(sid, ""),
                "category": "",
                "evidence_count": evc,
            })

    # Diseases per protein (dedupe by doid; rows can repeat per source)
    for row in conn.execute(
        "SELECT DISTINCT pd.uniprot_accession, d.doid, d.name, d.description "
        "FROM protein_disease pd "
        "JOIN disease d ON d.doid = pd.doid"
    ):
        ac = row["uniprot_accession"]
        if ac not in proteins_by_ac:
            continue
        bucket = proteins_by_ac[ac]["diseases"]
        if not any(e["disease_id"] == row["doid"] for e in bucket):
            bucket.append({
                "disease_id": row["doid"],
                "name": row["name"],
                "description": row["description"] or "",
            })

    conn.close()

    # Finalize each protein: counts + NLP augmentation
    empty_nlp = {f: (0 if f == "publication_count" else []) for f in NLP_FIELDS}
    proteins = []
    for ac, p in proteins_by_ac.items():
        p["glycosylation_count"] = len(p["glycosylation_sites"])
        p["phosphorylation_count"] = len(p["phosphorylation_sites"])
        p.update(aug.get(ac, empty_nlp))
        proteins.append(p)

    summary = {
        "total_proteins": len(proteins),
        "total_glyco_sites": sum(p["glycosylation_count"] for p in proteins),
        "total_phospho_sites": sum(p["phosphorylation_count"] for p in proteins),
        "total_diseases": len({d["disease_id"] for p in proteins for d in p["diseases"]}),
        "total_extracted_annotations": sum(len(p["ptm_annotations_extracted"]) for p in proteins),
        "total_snv_disease_links": sum(len(p["snv_with_disease_or_glycoeffect"]) for p in proteins),
        "data_source": "ptm_disease.db (live)",
        "cfde_integration": "GlyGen, MoTrPAC, GTEx, Kids First, IDG",
    }
    return {"summary": summary, "proteins": proteins}


def fetch_pqtl() -> dict:
    """Return {summary, overlaps[], protein_summaries{}} matching pqtl.json.

    Each overlap groups all PTM sites a single variant touches (direct or
    within 5 aa). Some fields that came from rich GlyGen SNV objects in the
    original build (chr_id, chr_pos, is_somatic/germline, variant_type,
    evidence_sources) aren't in the DB schema yet — those default to empty.
    """
    conn = get_conn()

    # Per-protein metadata
    protein_meta: dict[str, dict] = {}
    for row in conn.execute(
        "SELECT uniprot_accession, gene_symbol, protein_name FROM protein"
    ):
        ac = row["uniprot_accession"]
        gene = row["gene_symbol"] or ""
        pname = row["protein_name"] or ""
        display = f"{gene} ({pname})" if gene and pname else (gene or pname or ac)
        protein_meta[ac] = {"gene": gene, "name": display, "protein_name": pname}

    # Disease list per protein (overlap-level diseases inherit from protein).
    # protein_disease can have duplicate rows per source; dedupe by doid.
    diseases_by_protein: dict[str, list[dict]] = {}
    for row in conn.execute(
        "SELECT DISTINCT pd.uniprot_accession, d.doid, d.name "
        "FROM protein_disease pd JOIN disease d ON d.doid = pd.doid "
        "ORDER BY pd.uniprot_accession, d.doid"
    ):
        bucket = diseases_by_protein.setdefault(row["uniprot_accession"], [])
        if not any(e["id"] == row["doid"] for e in bucket):
            bucket.append({"id": row["doid"], "name": row["name"]})

    # All variant_ptm_site rows joined to variant + ptm_site + ptm_type
    rows = list(conn.execute("""
        SELECT v.variant_id, v.uniprot_accession, v.position AS snv_pos,
               v.ref_residue, v.alt_residue, v.dbsnp_id, v.clinical_significance,
               v.source_id AS variant_source,
               vps.distance, vps.effect, vps.notes,
               ps.position AS ptm_pos, ps.residue AS ptm_res,
               pt.name AS ptm_type_name
        FROM variant_ptm_site vps
        JOIN variant v ON v.variant_id = vps.variant_id
        JOIN ptm_site ps ON ps.ptm_site_id = vps.ptm_site_id
        JOIN ptm_type pt ON pt.ptm_type_id = ps.ptm_type_id
        ORDER BY v.uniprot_accession, v.position, ps.position
    """))

    # Per-protein totals for protein_summaries
    total_ptm_by_protein: dict[str, int] = {}
    for row in conn.execute(
        "SELECT uniprot_accession, COUNT(*) AS n FROM ptm_site GROUP BY uniprot_accession"
    ):
        total_ptm_by_protein[row["uniprot_accession"]] = row["n"]

    total_variants_by_protein: dict[str, int] = {}
    for row in conn.execute(
        "SELECT uniprot_accession, COUNT(*) AS n FROM variant GROUP BY uniprot_accession"
    ):
        total_variants_by_protein[row["uniprot_accession"]] = row["n"]

    conn.close()

    # Group by variant_id → one overlap entry per variant, with all PTM sites it touches
    by_variant: dict[int, dict] = {}
    for r in rows:
        vid = r["variant_id"]
        ac = r["uniprot_accession"]
        meta = protein_meta.get(ac, {"gene": "", "name": ac, "protein_name": ""})
        if vid not in by_variant:
            ref = r["ref_residue"] or ""
            by_variant[vid] = {
                "protein_ac": ac,
                "protein_name": meta["name"],
                "gene_name": meta["gene"],
                "snv_pos": r["snv_pos"],
                "snv_label": f"{ref}{r['snv_pos']}" if ref else str(r["snv_pos"]),
                "ref_aa": ref,
                "alt_aa": r["alt_residue"] or "",
                "chr_id": "",
                "chr_pos": "",
                "overlap_type": "direct" if r["distance"] == 0 else "proximal",
                "distance": r["distance"],
                "ptm_sites": [],
                "diseases": diseases_by_protein.get(ac, []),
                "glycoeffects": [],
                "is_disease_associated": bool(diseases_by_protein.get(ac)),
                "is_somatic": False,
                "is_germline": False,
                "variant_type": "",
                "comment": r["notes"] or "",
                "dbsnp": r["dbsnp_id"] or "",
                "evidence_sources": [r["variant_source"]] if r["variant_source"] else [],
            }
        entry = by_variant[vid]
        res1 = r["ptm_res"] or ""
        res3 = AA_3.get(res1, res1)
        entry["ptm_sites"].append({
            "pos": r["ptm_pos"],
            "label": f"{res3}{r['ptm_pos']}" if res3 else str(r["ptm_pos"]),
            "types": [r["ptm_type_name"].capitalize()] if r["ptm_type_name"] else [],
            "residue": res3,
        })
        if r["distance"] == 0:
            entry["overlap_type"] = "direct"
            entry["distance"] = 0
        else:
            if entry["overlap_type"] != "direct" and r["distance"] < entry["distance"]:
                entry["distance"] = r["distance"]
        if r["effect"] and r["effect"] != "unknown" and r["effect"] not in entry["glycoeffects"]:
            entry["glycoeffects"].append(r["effect"])

    overlaps = list(by_variant.values())

    # protein_summaries — only proteins that have at least one overlap
    protein_summaries: dict[str, dict] = {}
    for o in overlaps:
        ac = o["protein_ac"]
        ps = protein_summaries.setdefault(ac, {
            "name": o["protein_name"],
            "gene": o["gene_name"],
            "total_ptm_sites": total_ptm_by_protein.get(ac, 0),
            "total_snvs": total_variants_by_protein.get(ac, 0),
            "snvs_at_ptm": 0,
            "snvs_near_ptm": 0,
            "snvs_disease_and_ptm": 0,
        })
        if o["overlap_type"] == "direct":
            ps["snvs_at_ptm"] += 1
        else:
            ps["snvs_near_ptm"] += 1
        if o["is_disease_associated"]:
            ps["snvs_disease_and_ptm"] += 1

    summary = {
        "total_snv_ptm_overlaps": len(overlaps),
        "direct_overlaps": sum(1 for o in overlaps if o["overlap_type"] == "direct"),
        "proximal_overlaps": sum(1 for o in overlaps if o["overlap_type"] == "proximal"),
        "disease_associated": sum(1 for o in overlaps if o["is_disease_associated"]),
        "with_glycoeffect": sum(1 for o in overlaps if o["glycoeffects"]),
        "unique_diseases": len({d["id"] for o in overlaps for d in o["diseases"]}),
        "proteins_with_overlaps": len(protein_summaries),
    }
    return {"summary": summary, "overlaps": overlaps, "protein_summaries": protein_summaries}


def _ptm_type_label(ptm_name: str) -> tuple[str, str]:
    """Map ptm_type.name to (top-level type, subtype) for the site_evidence shape.

    'O-linked glycosylation' -> ('Glycosylation', 'O-linked')
    'phosphorylation'        -> ('Phosphorylation', '')
    """
    if not ptm_name:
        return ("", "")
    n = ptm_name.strip()
    if "glycosylation" in n.lower():
        subtype = n.lower().replace("glycosylation", "").strip()
        return ("Glycosylation", subtype.capitalize() if subtype else "")
    return (n.capitalize(), "")


def fetch_protein_evidence(uniprot_ac: str) -> dict:
    """Real-time PubMed-mined evidence for one protein.

    Mirrors the per-protein slice of site_evidence.json: PTM sites from the
    DB, PMIDs from the ptm_evidence table, abstracts fetched fresh from
    NCBI (cached locally), and disease/tissue/method regex mining applied
    to each abstract's title+text.
    """
    conn = get_conn()

    pinfo = conn.execute(
        "SELECT uniprot_accession, gene_symbol, protein_name "
        "FROM protein WHERE uniprot_accession=?",
        (uniprot_ac,),
    ).fetchone()
    if pinfo is None:
        conn.close()
        return {"error": f"protein {uniprot_ac} not in DB", "sites": []}

    gene = pinfo["gene_symbol"] or ""
    pname = pinfo["protein_name"] or ""
    display = f"{gene} ({pname})" if gene and pname else (gene or pname or uniprot_ac)

    # All PTM sites for this protein + their type info + glycan
    sites_raw = conn.execute("""
        SELECT ps.ptm_site_id, ps.position, ps.residue, ps.subtype AS ps_subtype,
               pt.name AS ptm_name, pt.category AS ptm_category
        FROM ptm_site ps
        JOIN ptm_type pt ON pt.ptm_type_id = ps.ptm_type_id
        WHERE ps.uniprot_accession = ?
        ORDER BY ps.position
    """, (uniprot_ac,)).fetchall()

    site_ids = [s["ptm_site_id"] for s in sites_raw]
    if not site_ids:
        conn.close()
        return {
            "uniprot_ac": uniprot_ac,
            "display_name": display,
            "gene_name": gene,
            "protein_name": pname,
            "sites": [],
            "summary": {"total_sites": 0, "unique_pmids": 0,
                        "abstracts_with_text": 0, "from_cache": 0},
        }

    placeholders = ",".join("?" * len(site_ids))
    pmids_by_site: dict[int, list[str]] = {sid: [] for sid in site_ids}
    seen_per_site: dict[int, set[str]] = {sid: set() for sid in site_ids}
    for row in conn.execute(
        f"SELECT ptm_site_id, pubmed_id FROM ptm_evidence "
        f"WHERE ptm_site_id IN ({placeholders}) AND pubmed_id IS NOT NULL",
        site_ids,
    ):
        sid, pmid = row["ptm_site_id"], str(row["pubmed_id"])
        if pmid not in seen_per_site[sid]:
            seen_per_site[sid].add(pmid)
            pmids_by_site[sid].append(pmid)

    glycan_by_site: dict[int, str] = {}
    for row in conn.execute(
        f"SELECT psm.ptm_site_id, m.modifier_id FROM ptm_site_modifier psm "
        f"JOIN modifier m ON m.modifier_id = psm.modifier_id "
        f"WHERE m.modifier_type='glycan' AND psm.ptm_site_id IN ({placeholders})",
        site_ids,
    ):
        glycan_by_site.setdefault(row["ptm_site_id"], row["modifier_id"])

    conn.close()

    # Fetch all unique PMIDs once (cache-first)
    all_pmids = sorted({p for plist in pmids_by_site.values() for p in plist})
    abstracts = pubmed.fetch_abstracts(all_pmids) if all_pmids else {}

    # Mine each abstract once
    per_pmid_annot: dict[str, dict] = {}
    abstracts_with_text = 0
    for pmid, rec in abstracts.items():
        title = rec.get("title", "")
        abstract = rec.get("abstract", "")
        if title or abstract:
            abstracts_with_text += 1
        diseases, tissues, methods = pubmed.mine_text(f"{title} {abstract}")
        per_pmid_annot[pmid] = {
            "title": title,
            "diseases": diseases,
            "tissues": tissues,
            "methods": methods,
        }

    # Assemble per-site records
    enriched: list[dict] = []
    for s in sites_raw:
        sid = s["ptm_site_id"]
        pos = s["position"]
        res1 = s["residue"] or ""
        res3 = AA_3.get(res1, res1)
        ptm_type, ptm_subtype = _ptm_type_label(s["ptm_name"])
        site_pmids = pmids_by_site.get(sid, [])

        site_diseases: set[str] = set()
        site_tissues: set[str] = set()
        site_methods: set[str] = set()
        evidence: list[dict] = []
        for pmid in site_pmids:
            ann = per_pmid_annot.get(pmid)
            if not ann:
                continue
            site_diseases.update(ann["diseases"])
            site_tissues.update(ann["tissues"])
            site_methods.update(ann["methods"])
            evidence.append({
                "p": pmid,
                "tt": ann["title"],
                "pmid": pmid,
                "title": ann["title"],
                "diseases": ann["diseases"],
                "tissues": ann["tissues"],
                "methods": ann["methods"],
            })

        enriched.append({
            "uniprot_ac": uniprot_ac,
            "display_name": display,
            "gene_name": gene,
            "protein_name": pname,
            "site_lbl": f"{res3}{pos}" if res3 else str(pos),
            "start_pos": pos,
            "residue": res3,
            "ptm_type": ptm_type,
            "ptm_subtype": ptm_subtype,
            "glycan_subtype": "",
            "glytoucan_ac": glycan_by_site.get(sid, ""),
            "category": "",
            "pmids": site_pmids,
            "conditions": sorted(site_diseases),
            "tissues": sorted(site_tissues),
            "methods": sorted(site_methods),
            "evidence": evidence,
        })

    return {
        "uniprot_ac": uniprot_ac,
        "display_name": display,
        "gene_name": gene,
        "protein_name": pname,
        "sites": enriched,
        "summary": {
            "total_sites": len(enriched),
            "unique_pmids": len(all_pmids),
            "abstracts_with_text": abstracts_with_text,
            "sites_with_condition_context": sum(1 for s in enriched if s["conditions"]),
            "sites_with_tissue_context": sum(1 for s in enriched if s["tissues"]),
        },
    }


def _protein_pmids(conn: sqlite3.Connection, uniprot_ac: str) -> list[str]:
    """All unique PubMed IDs evidenced for any PTM site on this protein."""
    rows = conn.execute("""
        SELECT DISTINCT e.pubmed_id
        FROM ptm_evidence e
        JOIN ptm_site ps ON ps.ptm_site_id = e.ptm_site_id
        WHERE ps.uniprot_accession = ? AND e.pubmed_id IS NOT NULL
    """, (uniprot_ac,)).fetchall()
    return [str(r[0]) for r in rows]


def fetch_protein_quantitative(uniprot_ac: str) -> dict:
    """Real-time sentence-level mining of quantitative PTM claims for one protein.

    Pulls all PMIDs evidenced for this protein's PTM sites, fetches abstracts
    (cached), splits each into sentences, and runs the mine_sentence regex
    extraction looking for site + direction + (disease|tissue|ptm).
    """
    conn = get_conn()
    pinfo = conn.execute(
        "SELECT uniprot_accession, gene_symbol, protein_name "
        "FROM protein WHERE uniprot_accession=?", (uniprot_ac,)
    ).fetchone()
    if pinfo is None:
        conn.close()
        return {"error": f"protein {uniprot_ac} not in DB", "hits": []}
    gene = pinfo["gene_symbol"] or ""
    pname = pinfo["protein_name"] or ""
    display = f"{gene} ({pname})" if gene and pname else (gene or pname or uniprot_ac)
    pmids = _protein_pmids(conn, uniprot_ac)
    conn.close()

    abstracts = pubmed.fetch_abstracts(pmids) if pmids else {}

    hits: list[dict] = []
    for pmid, rec in abstracts.items():
        title = rec.get("title", "")
        abstract = rec.get("abstract", "")
        full = f"{title}. {abstract}".strip()
        if not full or full == ".":
            continue
        for sent in pubmed.split_sentences(full):
            if len(sent) < 25:
                continue
            hit = pubmed.mine_sentence(sent, uniprot_ac, display, gene)
            if hit:
                hit["source"] = "PubMed abstract"
                hit["confidence"] = "medium"
                hit["pmid"] = pmid
                hit["paper_title"] = title[:150]
                hits.append(hit)

    summary = {
        "total_hits": len(hits),
        "unique_pmids_scanned": len(pmids),
        "abstracts_with_text": sum(1 for r in abstracts.values()
                                   if r.get("title") or r.get("abstract")),
        "with_direction": sum(1 for h in hits if h["direction"]),
        "with_quantitative": sum(1 for h in hits if h.get("quantitative")),
        "unique_conditions": len({d for h in hits for d in h["diseases"]}),
        "up_hits": sum(1 for h in hits if h["direction"] == "up"),
        "down_hits": sum(1 for h in hits if h["direction"] == "down"),
    }
    return {
        "uniprot_ac": uniprot_ac,
        "display_name": display,
        "gene_name": gene,
        "protein_name": pname,
        "hits": hits,
        "summary": summary,
    }


def _protein_snvs(conn: sqlite3.Connection, uniprot_ac: str) -> list[dict]:
    """All variants for this protein, each with its own BioMuta-style disease
    list (not the protein-wide bag), variant type, PTM overlaps, and comment.

    Output shape matches the per-protein 'snvs' field used by index.html's
    Track Viewer (compact keys: pos, lbl, ref, alt, ot, dist, ge, dis, vt,
    dbsnp, ptm, cmt).
    """
    # Per-variant disease annotations from GlyGen (cached to disk on first
    # call per protein). Falls back to {} if the API is unreachable —
    # variants then show no disease tags rather than wrong protein-level ones.
    glygen_by_pos = glygen.snv_disease_map(uniprot_ac)

    overlap_rows = list(conn.execute("""
        SELECT v.variant_id, v.position, v.ref_residue, v.alt_residue,
               v.dbsnp_id, v.clinical_significance,
               vps.distance, vps.effect,
               ps.position AS ptm_pos, ps.residue AS ptm_res
        FROM variant v
        LEFT JOIN variant_ptm_site vps ON vps.variant_id = v.variant_id
        LEFT JOIN ptm_site ps ON ps.ptm_site_id = vps.ptm_site_id
        WHERE v.uniprot_accession = ?
        ORDER BY v.position
    """, (uniprot_ac,)))

    by_variant: dict[int, dict] = {}
    for r in overlap_rows:
        vid = r["variant_id"]
        if vid not in by_variant:
            ref = r["ref_residue"] or ""
            pos = r["position"]
            ann = glygen_by_pos.get(pos, {})
            by_variant[vid] = {
                "pos": pos,
                "lbl": f"{ref}{pos}" if ref else str(pos),
                "ref": ref,
                "alt": r["alt_residue"] or "",
                "ot": None,
                "dist": None,
                "ge": [],
                "dis": ann.get("diseases", []),
                "vt": ann.get("variant_type", ""),
                "cmt": ann.get("comment", ""),
                "src": ann.get("sources", []),
                "chr_id": ann.get("chr_id", ""),
                "chr_pos": ann.get("chr_pos", ""),
                "dbsnp": r["dbsnp_id"] or "",
                "ptm": [],
            }
        if r["ptm_pos"] is not None:
            v = by_variant[vid]
            if v["ot"] is None or (r["distance"] is not None and r["distance"] == 0):
                v["ot"] = "direct" if r["distance"] == 0 else "proximal"
                v["dist"] = r["distance"]
            res3 = AA_3.get(r["ptm_res"] or "", r["ptm_res"] or "")
            v["ptm"].append({
                "pos": r["ptm_pos"],
                "lbl": f"{res3}{r['ptm_pos']}" if res3 else str(r["ptm_pos"]),
            })
            if r["effect"] and r["effect"] != "unknown" and r["effect"] not in v["ge"]:
                v["ge"].append(r["effect"])
    return list(by_variant.values())


def fetch_protein_tracks(uniprot_ac: str) -> dict:
    """Single-protein object matching the shape of tracks.json proteins[ac].

    Reuses fetch_protein_evidence + fetch_protein_quantitative — both share
    the PubMed cache so the second call is essentially free. SNVs come
    straight from the DB. Scoring follows build_scores.py:
      +1/-1 per hit, fold-change overrides, curated weight 2.
    """
    ev = fetch_protein_evidence(uniprot_ac)
    if ev.get("error"):
        return ev
    qt = fetch_protein_quantitative(uniprot_ac)

    # Index quantitative hits by site position AND aggregate at protein level
    quant_by_pos: dict[int, list[dict]] = {}
    for h in qt.get("hits", []):
        positions = []
        for s in h.get("sites", []):
            m = re.search(r"\d+", s)
            if m:
                positions.append(int(m.group(0)))
        for p in positions:
            quant_by_pos.setdefault(p, []).append(h)

    # Helper: per-hit score (matches build_scores.py)
    def hit_score(h: dict) -> float:
        d = h.get("direction")
        quant = h.get("quantitative")
        conf = h.get("confidence", "medium")
        if quant and quant.get("value") is not None:
            val = quant["value"] / 100 if quant.get("type") == "percent" else quant["value"]
            return val if d == "up" else -val
        base = 2.0 if conf == "high" else 1.0
        return base if d == "up" else -base

    # Protein-level cumulative score per condition — uses every hit (with
    # or without site numbers), so sites that have a condition mentioned in
    # their evidence but no site-specific direction get a non-zero fallback.
    protein_cond_scores: dict[str, dict] = {}
    for h in qt.get("hits", []):
        s = hit_score(h)
        for cond in h.get("diseases", []):
            if cond in ("Control", "Normal/healthy"):
                continue
            pcs = protein_cond_scores.setdefault(cond, {
                "score": 0.0, "up": 0, "down": 0, "evidence": [],
            })
            if h.get("direction") == "up":
                pcs["up"] += 1
            else:
                pcs["down"] += 1
            pcs["score"] += s
            if len(pcs["evidence"]) < 3:
                pcs["evidence"].append({
                    "pmid": h.get("pmid", ""),
                    "sentence": h.get("sentence", "")[:150],
                    "score": s,
                })

    # Build sites list in the tracks shape
    sites_out: list[dict] = []
    all_conditions: set[str] = set()
    n_sites_total = max(len(ev.get("sites", [])), 1)

    for ev_site in ev.get("sites", []):
        pos = ev_site["start_pos"]
        res3 = ev_site["residue"]
        ptm_type = ev_site["ptm_type"]
        ptm_subtype = ev_site["ptm_subtype"]
        ptm_types = []
        if ptm_type:
            label = f"{ptm_type} ({ptm_subtype})" if ptm_subtype else ptm_type
            ptm_types.append(label)
        glytoucan = [ev_site["glytoucan_ac"]] if ev_site.get("glytoucan_ac") else []

        # Conditions from evidence (observed-only, count = #PMIDs mentioning)
        cond_dict: dict[str, dict] = {}
        for ev_rec in ev_site.get("evidence", []):
            for c in ev_rec.get("diseases", []):
                e = cond_dict.setdefault(c, {"direction": "observed", "count": 0, "evidence": []})
                e["count"] += 1

        # Scores + direction overrides from quantitative hits at this pos
        scores: dict[str, dict] = {}
        for h in quant_by_pos.get(pos, []):
            s = hit_score(h)
            for cond in h.get("diseases", []):
                if cond in ("Control", "Normal/healthy"):
                    continue
                scs = scores.setdefault(cond, {
                    "score": 0.0, "up": 0, "down": 0, "max_fc": None, "evidence": [],
                })
                if h.get("direction") == "up":
                    scs["up"] += 1
                else:
                    scs["down"] += 1
                scs["score"] += s
                q = h.get("quantitative")
                if q and q.get("value") is not None:
                    cand = q["value"] if h.get("direction") == "up" else -q["value"]
                    if scs["max_fc"] is None or abs(cand) > abs(scs["max_fc"]):
                        scs["max_fc"] = cand
                if len(scs["evidence"]) < 3:
                    scs["evidence"].append({
                        "pmid": h.get("pmid", ""),
                        "sentence": h.get("sentence", "")[:150],
                        "score": s,
                    })
                # Reflect into cond_dict direction
                ce = cond_dict.setdefault(cond, {"direction": "observed", "count": 0, "evidence": []})
                ce["direction"] = h.get("direction") or ce["direction"]
                if len(ce["evidence"]) < 2:
                    ce["evidence"].append({
                        "sentence": h.get("sentence", "")[:150],
                        "pmid": h.get("pmid", ""),
                        "confidence": h.get("confidence", "medium"),
                        "quantitative": h.get("quantitative"),
                    })

        # Protein-level fallback: for every condition mentioned in this
        # site's evidence that doesn't already have a site-level score, use
        # the diluted protein-level total (matches build_scores.py logic).
        for cond in cond_dict:
            if cond in scores or cond not in protein_cond_scores:
                continue
            pcs = protein_cond_scores[cond]
            scores[cond] = {
                "score": pcs["score"] / n_sites_total,
                "up": pcs["up"],
                "down": pcs["down"],
                "max_fc": None,
                "evidence": list(pcs["evidence"]),
                "protein_level": True,
            }

        # Round scores for display
        for c, v in scores.items():
            v["score"] = round(v["score"], 2)
            if v["max_fc"] is not None:
                v["max_fc"] = round(v["max_fc"], 2)

        all_conditions.update(cond_dict.keys())
        all_conditions.update(scores.keys())

        sites_out.append({
            "pos": pos,
            "residue": res3,
            "label": ev_site["site_lbl"],
            "ptm_types": ptm_types,
            "glytoucan": glytoucan[:3],
            "tissues": ev_site.get("tissues", []),
            "conditions": cond_dict,
            "scores": scores,
        })

    # SNVs from DB
    conn = get_conn()
    snvs = _protein_snvs(conn, uniprot_ac)
    pinfo = conn.execute(
        "SELECT sequence_length FROM protein WHERE uniprot_accession=?",
        (uniprot_ac,),
    ).fetchone()
    length = pinfo["sequence_length"] if pinfo and pinfo["sequence_length"] else 0
    diseases_rows = conn.execute(
        "SELECT DISTINCT d.name FROM protein_disease pd "
        "JOIN disease d ON d.doid = pd.doid WHERE pd.uniprot_accession=? "
        "ORDER BY d.name",
        (uniprot_ac,),
    ).fetchall()
    conn.close()
    disease_names = [r["name"] for r in diseases_rows]

    return {
        "ac": uniprot_ac,
        "name": ev["display_name"],
        "gene": ev["gene_name"],
        "protein_name": ev["protein_name"],
        "length": length,
        "seq_snippet": "",
        "sites": sites_out,
        "total_sites": len(sites_out),
        "conditions": sorted(all_conditions),
        "diseases": disease_names,
        "snvs": snvs,
    }


def fetch_tracks_manifest() -> dict:
    """Lightweight protein manifest for the Track Viewer dropdown.

    Returns proteins keyed by ac with just ac/name/gene/total_sites/length —
    no sites, no scores. The actual per-protein track is loaded lazily on
    selection via fetch_protein_tracks(ac).
    """
    conn = get_conn()
    counts = {r["uniprot_accession"]: r["n"] for r in conn.execute(
        "SELECT uniprot_accession, COUNT(*) AS n FROM ptm_site GROUP BY uniprot_accession"
    )}
    proteins: dict[str, dict] = {}
    for row in conn.execute(
        "SELECT uniprot_accession, gene_symbol, protein_name, sequence_length "
        "FROM protein ORDER BY gene_symbol, uniprot_accession"
    ):
        ac = row["uniprot_accession"]
        n = counts.get(ac, 0)
        if n == 0:
            continue
        gene = row["gene_symbol"] or ""
        pname = row["protein_name"] or ""
        display = f"{gene} ({pname})" if gene and pname else (gene or pname or ac)
        proteins[ac] = {
            "ac": ac,
            "name": display,
            "gene": gene,
            "protein_name": pname,
            "length": row["sequence_length"] or 0,
            "total_sites": n,
            # sites/conditions/scores/snvs are lazy — fetched on demand
            "sites": [],
            "conditions": [],
            "diseases": [],
            "snvs": [],
            "_lazy": True,
        }
    conn.close()
    return {
        "conditions": [],  # populated per-protein after lazy fetch
        "proteins": proteins,
        "scoring_method": {
            "description": "Evidence-weighted literature score per site per condition (live)",
            "up_abstract": "+1 per PubMed abstract reporting increase",
            "down_abstract": "-1 per PubMed abstract reporting decrease",
            "up_curated": "+2 per curated UniProtKB annotation reporting increase",
            "down_curated": "-2 per curated UniProtKB annotation reporting decrease",
            "fold_change": "If abstract mentions explicit fold-change, that value is used instead",
            "interpretation": "Positive score = evidence of PTM increase; negative = decrease; magnitude = evidence strength",
        },
    }
