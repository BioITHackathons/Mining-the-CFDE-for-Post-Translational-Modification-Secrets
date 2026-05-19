#!/usr/bin/env python3
"""
Main orchestrator: initialise the PTM-disease database, seed reference data,
then run all loaders.

Usage:
    python loaders/load_all.py                 # full pipeline
    python loaders/load_all.py --skip-glygen   # skip GlyGen direct API (slow)
    python loaders/load_all.py --db /path/to.db  # custom DB path
"""

import argparse
import sys
from pathlib import Path

# Ensure loaders/ is on the path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from db_utils import init_db, get_connection, DB_PATH, upsert_data_source, upsert_protein, upsert_disease, upsert_protein_disease
from seed_proteins import SEED_PROTEINS, SEED_DISEASES


def seed_reference_data(conn):
    """Insert seed proteins, diseases, and their associations."""
    print("Seeding reference proteins and diseases...")

    # Proteins
    for acc, gene in SEED_PROTEINS.items():
        upsert_protein(conn, acc, gene)

    # Diseases and protein-disease links
    for doid, info in SEED_DISEASES.items():
        upsert_disease(conn, doid, info["name"], category=info["category"])
        for acc in info["proteins"]:
            upsert_protein_disease(conn, acc, doid, association_type="genetic", source_id=None)

    conn.commit()
    print(f"  Seeded {len(SEED_PROTEINS)} proteins, {len(SEED_DISEASES)} diseases.")


def run_glygen_loader(conn):
    """Run the GlyGen direct API loader."""
    print("\n" + "=" * 60)
    print("LOADER: GlyGen Direct API")
    print("=" * 60)
    import load_glygen
    accessions = list(SEED_PROTEINS.keys())
    load_glygen.main(accessions=accessions)


def run_cfde_kc_loader(conn):
    """Run the CFDE Knowledge Center loader."""
    print("\n" + "=" * 60)
    print("LOADER: CFDE Knowledge Center (GTEx, Kids First, IDG, GlyGen KC)")
    print("=" * 60)
    import load_cfde_kc
    genes = list({v for v in SEED_PROTEINS.values()})
    load_cfde_kc.main(genes=genes)


def print_summary(conn):
    """Print database summary statistics."""
    print("\n" + "=" * 60)
    print("DATABASE SUMMARY")
    print("=" * 60)

    tables = [
        ("protein", "Proteins"),
        ("disease", "Diseases"),
        ("ptm_type", "PTM types"),
        ("data_source", "Data sources"),
        ("ptm_site", "PTM sites"),
        ("ptm_evidence", "Evidence records"),
        ("protein_disease", "Protein-disease links"),
        ("ptm_site_disease", "Site-disease links"),
        ("modifier", "Modifiers (glycans, etc.)"),
        ("ptm_site_modifier", "Site-modifier links"),
        ("variant", "Variants"),
        ("variant_ptm_site", "Variant-PTM site links"),
        ("tissue_expression", "Tissue expression records"),
        ("biomarker", "Biomarkers"),
        ("ptm_crosstalk", "PTM crosstalk records"),
    ]

    for table, label in tables:
        try:
            count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            print(f"  {label:.<40} {count:>8,}")
        except Exception:
            print(f"  {label:.<40} {'N/A':>8}")

    # PTM sites by type
    print("\n  PTM sites by type:")
    rows = conn.execute("""
        SELECT pt.name, COUNT(*) as n
        FROM ptm_site ps JOIN ptm_type pt ON pt.ptm_type_id = ps.ptm_type_id
        GROUP BY pt.name ORDER BY n DESC
    """).fetchall()
    for row in rows:
        print(f"    {row[0]:.<38} {row[1]:>8,}")

    # Proteins by disease
    print("\n  Proteins with PTM data by disease:")
    rows = conn.execute("""
        SELECT d.name, COUNT(DISTINCT ps.uniprot_accession) as n
        FROM ptm_site ps
        JOIN protein_disease pd ON pd.uniprot_accession = ps.uniprot_accession
        JOIN disease d ON d.doid = pd.doid
        GROUP BY d.name ORDER BY n DESC
    """).fetchall()
    for row in rows:
        print(f"    {row[0]:.<38} {row[1]:>8,}")


def main():
    parser = argparse.ArgumentParser(description="Build the PTM-disease database")
    parser.add_argument("--db", type=str, default=str(DB_PATH), help="Database file path")
    parser.add_argument("--skip-glygen", action="store_true", help="Skip GlyGen direct API loader")
    parser.add_argument("--skip-kc", action="store_true", help="Skip CFDE KC loader")
    args = parser.parse_args()

    db_path = Path(args.db)
    if db_path.exists():
        db_path.unlink()
        print(f"Removed existing database: {db_path}")

    print(f"Initialising database at {db_path}...")
    conn = init_db(db_path)

    # Seed
    seed_reference_data(conn)

    # Run loaders
    if not args.skip_glygen:
        run_glygen_loader(conn)

    if not args.skip_kc:
        run_cfde_kc_loader(conn)

    # Summary
    conn = get_connection(db_path)
    print_summary(conn)
    conn.close()

    print(f"\nDone. Database: {db_path}")
    print(f"Size: {db_path.stat().st_size / 1024:.1f} KB")


if __name__ == "__main__":
    main()
