"""
Loader: GlyGen direct API (https://api.glygen.org)
Extracts glycosylation, phosphorylation, SNVs, and disease associations
for a list of proteins.

Usage:
    python loaders/load_glygen.py            # uses default gene list
    python loaders/load_glygen.py P05067 P10636  # specific UniProt accessions
"""

import sys
import time
import json
import requests
from db_utils import (
    get_connection, upsert_data_source, upsert_protein, upsert_disease,
    upsert_ptm_type, upsert_ptm_site, insert_evidence, upsert_protein_disease,
    upsert_modifier, upsert_ptm_site_modifier, upsert_variant, TODAY,
)

API_BASE = "https://api.glygen.org"
SOURCE_ID = "glygen"

# ── PTM type mapping ─────────────────────────────────────────────────────────

GLYCO_TYPE_MAP = {
    "N-linked": ("MOD:00006", "N-linked glycosylation", "glycosylation"),
    "O-linked": ("MOD:00007", "O-linked glycosylation", "glycosylation"),
    "C-linked": ("MOD:01084", "C-linked glycosylation", "glycosylation"),
}
PHOSPHO_TYPE_ID = "MOD:00696"

# Three-letter to one-letter amino acid codes
AA_MAP = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def fetch_protein_detail(accession: str) -> dict | None:
    """Call GlyGen protein detail API."""
    url = f"{API_BASE}/protein/detail/{accession}"
    try:
        r = requests.post(url, json={}, timeout=60)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        print(f"  [WARN] Failed to fetch {accession}: {e}")
        return None


def _extract_pubmed_ids(evidence_list: list) -> list[str]:
    """Pull PubMed IDs from GlyGen evidence array."""
    return [e["id"] for e in evidence_list if e.get("database") == "PubMed"]


def _evidence_type_from_evidence(evidence_list: list) -> str:
    """Infer evidence type from GlyGen evidence sources."""
    dbs = {e.get("database", "").lower() for e in evidence_list}
    if "pubmed" in dbs:
        return "experimental"
    if "uniprotkb" in dbs or "iptmnet" in dbs:
        return "curated"
    return "inferred"


def load_protein(conn, accession: str):
    """Load all PTM data for one protein from GlyGen."""
    print(f"  Fetching {accession} from GlyGen...")
    data = fetch_protein_detail(accession)
    if not data:
        return

    # ── Protein record ────────────────────────────────────────────────────
    gene_names = data.get("gene_names", [])
    gene_symbol = gene_names[0].get("name", accession) if gene_names else accession
    protein_names = data.get("protein_names", [])
    protein_name = protein_names[0].get("name") if protein_names else None
    seq_len = len(data.get("sequence", {}).get("sequence", "")) or None

    upsert_protein(conn, accession, gene_symbol, protein_name, sequence_length=seq_len)

    # ── Diseases ──────────────────────────────────────────────────────────
    for d in data.get("disease", []):
        doid = d.get("disease_id", "")
        rn = d.get("recommended_name", {})
        if doid.startswith("DOID:"):
            upsert_disease(conn, doid, rn.get("name", doid), rn.get("description"))
            upsert_protein_disease(conn, accession, doid, source_id=SOURCE_ID)

    # ── Glycosylation sites ───────────────────────────────────────────────
    glyco_count = 0
    for site in data.get("glycosylation", []):
        gtype = site.get("type", "O-linked")
        ptm_type_id, ptm_name, ptm_cat = GLYCO_TYPE_MAP.get(
            gtype, ("MOD:00007", gtype + " glycosylation", "glycosylation")
        )
        upsert_ptm_type(conn, ptm_type_id, ptm_name, ptm_cat)

        pos = site.get("start_pos")
        residue_3 = site.get("residue") or site.get("start_aa")
        residue = AA_MAP.get(residue_3, residue_3)
        if not pos:
            continue

        site_id = upsert_ptm_site(
            conn, accession, pos, residue, ptm_type_id,
            subtype=site.get("subtype"),
            flanking_sequence=site.get("site_seq"),
        )

        # Evidence
        ev_list = site.get("evidence", [])
        pmids = _extract_pubmed_ids(ev_list)
        ev_type = _evidence_type_from_evidence(ev_list)
        for pmid in pmids:
            insert_evidence(conn, site_id, SOURCE_ID, ev_type, pubmed_id=pmid, raw_record=site)
        if not pmids:
            insert_evidence(conn, site_id, SOURCE_ID, ev_type, raw_record=site)

        # Glycan modifier
        glytoucan = site.get("glytoucan_ac")
        if glytoucan:
            upsert_modifier(conn, glytoucan, "glycan")
            upsert_ptm_site_modifier(
                conn, site_id, glytoucan, source_id=SOURCE_ID,
                evidence_type=ev_type,
                pubmed_id=pmids[0] if pmids else None,
            )
        glyco_count += 1

    # ── Phosphorylation sites ─────────────────────────────────────────────
    phos_count = 0
    upsert_ptm_type(conn, PHOSPHO_TYPE_ID, "phosphorylation", "phosphorylation")
    for site in data.get("phosphorylation", []):
        pos = site.get("start_pos")
        residue_3 = site.get("residue")
        residue = AA_MAP.get(residue_3, residue_3)
        if not pos:
            continue

        site_id = upsert_ptm_site(
            conn, accession, pos, residue, PHOSPHO_TYPE_ID,
            notes=site.get("comment"),
        )

        ev_list = site.get("evidence", [])
        pmids = _extract_pubmed_ids(ev_list)
        ev_type = _evidence_type_from_evidence(ev_list)
        for pmid in pmids:
            insert_evidence(conn, site_id, SOURCE_ID, ev_type, pubmed_id=pmid, raw_record=site)
        if not pmids:
            insert_evidence(conn, site_id, SOURCE_ID, ev_type, raw_record=site)
        phos_count += 1

    # ── SNVs ──────────────────────────────────────────────────────────────
    snv_count = 0
    for snv in data.get("snv", []):
        pos = snv.get("start_pos")
        if not pos:
            continue

        # Extract dbSNP ID
        dbsnp = None
        for ev in snv.get("evidence", []):
            if ev.get("database") == "dbSNP":
                dbsnp = ev.get("id")
                break

        variant_id = upsert_variant(
            conn, accession, pos,
            ref_residue=snv.get("sequence_org"),
            alt_residue=snv.get("sequence_mut"),
            dbsnp_id=dbsnp,
            source_id=SOURCE_ID,
        )

        # Link variant to any PTM sites at that position
        rows = conn.execute(
            "SELECT ptm_site_id FROM ptm_site WHERE uniprot_accession = ? AND position = ?",
            (accession, pos),
        ).fetchall()
        for row in rows:
            glycoeffects = snv.get("glycoeffect", [])
            effect = glycoeffects[0] if glycoeffects else "unknown"
            conn.execute(
                """INSERT OR IGNORE INTO variant_ptm_site (variant_id, ptm_site_id, effect, distance, notes)
                   VALUES (?, ?, ?, 0, ?)""",
                (variant_id, row[0], effect, snv.get("comment")),
            )

        # Link variant to diseases mentioned in the SNV record
        for d in snv.get("disease", []):
            doid = d.get("disease_id", "")
            if doid.startswith("DOID:"):
                rn = d.get("recommended_name", {})
                upsert_disease(conn, doid, rn.get("name", doid), rn.get("description"))
        snv_count += 1

    print(f"    -> {glyco_count} glyco, {phos_count} phospho, {snv_count} SNVs")


def main(accessions: list[str] | None = None):
    conn = get_connection()

    # Register GlyGen as a data source
    upsert_data_source(
        conn, SOURCE_ID, "GlyGen",
        url="https://glygen.org",
        cfde_program="GlyGen",
        description="Glycosylation and PTM data from GlyGen",
    )

    # Determine which accessions to load
    if accessions is None:
        if len(sys.argv) > 1:
            accessions = sys.argv[1:]
        else:
            from seed_proteins import SEED_PROTEINS
            accessions = list(SEED_PROTEINS.keys())

    print(f"Loading {len(accessions)} proteins from GlyGen...")
    for i, acc in enumerate(accessions):
        load_protein(conn, acc)
        conn.commit()
        if i < len(accessions) - 1:
            time.sleep(0.5)  # rate limit

    print(f"GlyGen load complete.")
    conn.close()


if __name__ == "__main__":
    main()
