"""
Loader: MoTrPAC (Molecular Transducers of Physical Activity Consortium)
Downloads training-regulated PTM features (phospho, acetyl, ubiquityl)
from the MotrpacRatTraining6moData R package (GitHub), maps rat sites
to human UniProt sites, and loads into the PTM-disease database.

Data source: https://github.com/MoTrPAC/MotrpacRatTraining6moData
Publication: Nature 629, 174-183 (2024). doi:10.1038/s41586-023-06877-w

Usage:
    python loaders/load_motrpac.py
"""

import sys
import re
import tempfile
import os
import requests
import pyreadr
import sqlite3
from pathlib import Path

# Allow importing from the loaders directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from db_utils import (
    get_connection, upsert_data_source, upsert_protein, upsert_ptm_type,
    upsert_ptm_site, insert_evidence, TODAY,
)
from seed_proteins import SEED_PROTEINS

SOURCE_ID = "motrpac"
DATA_BASE = "https://github.com/MoTrPAC/MotrpacRatTraining6moData/raw/main/data"

# PTM assay -> (ptm_type_id, ptm_name, category)
ASSAY_PTM_MAP = {
    "PHOSPHO": ("MOD:00696", "phosphorylation", "phosphorylation"),
    "ACETYL":  ("MOD:00394", "acetylation", "acetylation"),
    "UBIQ":    ("MOD:01148", "ubiquitylation", "ubiquitylation"),
}

# Parse site info from feature_ID like "NP_001003673.1_K477k"
# Format: RefSeqAcc_<Residue><Position><residue_lower>
SITE_RE = re.compile(r"^(.+?)_([A-Z])(\d+)([a-z])(.*)$")


def download_rda(filename: str) -> dict:
    """Download an .rda file from MoTrPAC GitHub and return as dict of DataFrames."""
    url = f"{DATA_BASE}/{filename}"
    print(f"  Downloading {filename}...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()

    tmp = tempfile.NamedTemporaryFile(suffix=".rda", delete=False)
    tmp.write(r.content)
    tmp.close()

    try:
        result = pyreadr.read_r(tmp.name)
        return result
    finally:
        os.unlink(tmp.name)


def parse_ptm_site(feature_id: str):
    """Parse a MoTrPAC feature_ID into (refseq_acc, residue, position).
    
    Examples:
        NP_001003673.1_K477k -> (NP_001003673.1, K, 477)
        NP_036620.2_S234sS237s -> (NP_036620.2, S, 234)  # first site only
    """
    m = SITE_RE.match(feature_id)
    if m:
        return m.group(1), m.group(2), int(m.group(3))
    return None, None, None


def main(genes: list[str] | None = None):
    conn = get_connection()

    # Register source
    upsert_data_source(
        conn, SOURCE_ID, "MoTrPAC",
        url="https://motrpac-data.org",
        cfde_program="MoTrPAC",
        description="PTM data from MoTrPAC endurance training study in rats, mapped to human orthologs. Nature 629, 174-183 (2024).",
    )

    # Register PTM types
    for assay, (ptm_id, ptm_name, ptm_cat) in ASSAY_PTM_MAP.items():
        upsert_ptm_type(conn, ptm_id, ptm_name, ptm_cat)
    conn.commit()

    # ── Download MoTrPAC data ─────────────────────────────────────────────
    print("Downloading MoTrPAC data files...")
    
    tr_features_data = download_rda("TRAINING_REGULATED_FEATURES.rda")
    tr_features = tr_features_data["TRAINING_REGULATED_FEATURES"]

    f2g_data = download_rda("FEATURE_TO_GENE_FILT.rda")
    f2g = f2g_data["FEATURE_TO_GENE_FILT"]

    r2h_gene_data = download_rda("RAT_TO_HUMAN_GENE.rda")
    r2h_gene = r2h_gene_data["RAT_TO_HUMAN_GENE"]

    r2h_phospho_data = download_rda("RAT_TO_HUMAN_PHOSPHO.rda")
    r2h_phospho = r2h_phospho_data["RAT_TO_HUMAN_PHOSPHO"]

    # ── Build lookup maps ─────────────────────────────────────────────────
    print("Building lookup maps...")

    # Rat gene symbol -> human gene symbol
    rat_to_human_gene = {}
    for _, row in r2h_gene.iterrows():
        rat_sym = row.get("RAT_SYMBOL")
        human_sym = row.get("HUMAN_ORTHOLOG_SYMBOL")
        if rat_sym and human_sym and str(human_sym) != "nan":
            rat_to_human_gene[str(rat_sym).upper()] = str(human_sym).upper()

    # Rat RefSeq PTM ID -> human UniProt PTM ID
    rat_to_human_ptm = {}
    for _, row in r2h_phospho.iterrows():
        rat_id = row.get("ptm_id_rat_refseq")
        human_id = row.get("ptm_id_human_uniprot")
        if rat_id and human_id and str(human_id) != "nan":
            rat_to_human_ptm[str(rat_id)] = str(human_id)

    # Feature ID -> rat gene symbol
    feature_to_gene = {}
    for _, row in f2g.iterrows():
        fid = row.get("feature_ID")
        gene = row.get("gene_symbol")
        if fid and gene:
            feature_to_gene[str(fid)] = str(gene)

    # Human gene symbol -> UniProt accession (from our seed proteins)
    human_gene_to_acc = {g: a for a, g in SEED_PROTEINS.items()}

    # Determine target genes
    if genes is None:
        target_genes = set(SEED_PROTEINS.values())
    else:
        target_genes = set(g.upper() for g in genes)

    # ── Filter to PTM assays ──────────────────────────────────────────────
    ptm_features = tr_features[tr_features["assay"].isin(["PHOSPHO", "ACETYL", "UBIQ"])].copy()
    print(f"Total training-regulated PTM features: {len(ptm_features)}")
    print(f"  PHOSPHO: {len(ptm_features[ptm_features['assay']=='PHOSPHO'])}")
    print(f"  ACETYL:  {len(ptm_features[ptm_features['assay']=='ACETYL'])}")
    print(f"  UBIQ:    {len(ptm_features[ptm_features['assay']=='UBIQ'])}")

    # ── Process features ──────────────────────────────────────────────────
    # Deduplicate to unique (feature_ID, assay, tissue) -- aggregate across timepoints later
    unique_sites = ptm_features.drop_duplicates(subset=["feature_ID", "assay", "tissue"])
    print(f"Unique (feature, assay, tissue) combinations: {len(unique_sites)}")

    loaded = 0
    skipped_no_gene = 0
    skipped_no_human = 0
    skipped_no_acc = 0
    skipped_no_site = 0

    seen_sites = set()

    for _, row in unique_sites.iterrows():
        fid = str(row["feature_ID"])
        assay = row["assay"]
        tissue = row["tissue"]

        # Get rat gene symbol
        rat_gene = feature_to_gene.get(fid)
        if not rat_gene:
            skipped_no_gene += 1
            continue

        # Map to human gene
        human_gene = rat_to_human_gene.get(rat_gene.upper())
        if not human_gene:
            skipped_no_human += 1
            continue

        # Check if this gene is in our target set
        if human_gene.upper() not in target_genes:
            continue

        # Get UniProt accession
        acc = human_gene_to_acc.get(human_gene.upper())
        if not acc:
            skipped_no_acc += 1
            continue

        # Parse site from feature_ID
        ptm_type_id, ptm_name, ptm_cat = ASSAY_PTM_MAP[assay]

        # For phospho, try rat-to-human PTM site mapping
        position = None
        residue = None

        if assay == "PHOSPHO":
            human_ptm_id = rat_to_human_ptm.get(fid)
            if human_ptm_id:
                # Parse human UniProt PTM ID like "P05067_S198s"
                parts = human_ptm_id.split("_", 1)
                if len(parts) == 2:
                    m = re.match(r"([A-Z])(\d+)", parts[1])
                    if m:
                        residue = m.group(1)
                        position = int(m.group(2))

        # Fallback: parse from rat feature_ID (position may differ from human)
        if position is None:
            _, rat_residue, rat_pos = parse_ptm_site(fid)
            if rat_pos:
                residue = rat_residue
                position = rat_pos
            else:
                skipped_no_site += 1
                continue

        # Deduplicate: only insert each (acc, position, ptm_type) once
        site_key = (acc, position, ptm_type_id)
        if site_key in seen_sites:
            # Still add evidence for this tissue
            continue
        seen_sites.add(site_key)

        upsert_protein(conn, acc, human_gene)

        site_id = upsert_ptm_site(
            conn, acc, position, residue, ptm_type_id,
            notes=f"MoTrPAC training-regulated; rat gene={rat_gene}",
        )

        # Get summary stats across all timepoints for this feature/tissue
        feature_rows = ptm_features[
            (ptm_features["feature_ID"] == fid) &
            (ptm_features["tissue"] == tissue)
        ]
        min_p = feature_rows["timewise_p_value"].min()
        max_lfc = feature_rows["timewise_logFC"].abs().max()
        training_q = feature_rows["training_q"].iloc[0] if len(feature_rows) > 0 else None

        insert_evidence(
            conn, site_id, SOURCE_ID, "experimental",
            method="mass spectrometry (TMT)",
            confidence_score=1 - min_p if min_p is not None else None,
            raw_record={
                "feature_ID": fid,
                "assay": assay,
                "tissue": tissue,
                "rat_gene": rat_gene,
                "human_gene": human_gene,
                "min_timewise_p": float(min_p) if min_p is not None else None,
                "max_abs_logFC": float(max_lfc) if max_lfc is not None else None,
                "training_q": float(training_q) if training_q is not None else None,
                "organism": "Rattus norvegicus (mapped to human)",
            },
        )
        loaded += 1

    conn.commit()

    print(f"\nMoTrPAC load complete:")
    print(f"  PTM sites loaded: {loaded}")
    print(f"  Skipped (no gene mapping): {skipped_no_gene}")
    print(f"  Skipped (no human ortholog): {skipped_no_human}")
    print(f"  Skipped (no UniProt acc): {skipped_no_acc}")
    print(f"  Skipped (no site parsed): {skipped_no_site}")

    conn.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(genes=sys.argv[1:])
    else:
        main()
