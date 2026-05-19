"""
Shared database utilities for PTM-disease database loaders.
"""

import sqlite3
import json
from pathlib import Path
from datetime import date

DB_PATH = Path(__file__).resolve().parent.parent / "ptm_disease.db"
SCHEMA_PATH = Path(__file__).resolve().parent.parent / "ptm_disease_db_schema.sql"

TODAY = date.today().isoformat()


def get_connection(db_path: str | Path | None = None) -> sqlite3.Connection:
    """Return a connection with WAL mode and foreign keys enabled."""
    path = str(db_path or DB_PATH)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA foreign_keys = ON")
    conn.row_factory = sqlite3.Row
    return conn


def init_db(db_path: str | Path | None = None) -> sqlite3.Connection:
    """Create the database from the schema DDL. Returns the connection."""
    conn = get_connection(db_path)
    schema_sql = SCHEMA_PATH.read_text()
    conn.executescript(schema_sql)
    conn.commit()
    return conn


# ── Upsert helpers ───────────────────────────────────────────────────────────


def upsert_data_source(conn, source_id, name, url=None, cfde_program=None, description=None):
    conn.execute(
        """INSERT OR IGNORE INTO data_source (source_id, name, url, cfde_program, description)
           VALUES (?, ?, ?, ?, ?)""",
        (source_id, name, url, cfde_program, description),
    )


def upsert_protein(conn, accession, gene_symbol, protein_name=None, organism="Homo sapiens",
                    sequence_length=None, reviewed=1):
    conn.execute(
        """INSERT INTO protein (uniprot_accession, gene_symbol, protein_name, organism, sequence_length, reviewed)
           VALUES (?, ?, ?, ?, ?, ?)
           ON CONFLICT(uniprot_accession) DO UPDATE SET
             gene_symbol = COALESCE(excluded.gene_symbol, protein.gene_symbol),
             protein_name = COALESCE(excluded.protein_name, protein.protein_name),
             sequence_length = COALESCE(excluded.sequence_length, protein.sequence_length)""",
        (accession, gene_symbol, protein_name, organism, sequence_length, reviewed),
    )


def upsert_disease(conn, doid, name, description=None, category=None):
    conn.execute(
        """INSERT INTO disease (doid, name, description, category)
           VALUES (?, ?, ?, ?)
           ON CONFLICT(doid) DO UPDATE SET
             name = COALESCE(excluded.name, disease.name),
             description = COALESCE(excluded.description, disease.description),
             category = COALESCE(excluded.category, disease.category)""",
        (doid, name, description, category),
    )


def upsert_ptm_type(conn, ptm_type_id, name, category=None, description=None):
    conn.execute(
        """INSERT OR IGNORE INTO ptm_type (ptm_type_id, name, category, description)
           VALUES (?, ?, ?, ?)""",
        (ptm_type_id, name, category, description),
    )


def upsert_ptm_site(conn, accession, position, residue, ptm_type_id,
                     subtype=None, flanking_sequence=None, notes=None):
    """Insert a PTM site; returns the ptm_site_id."""
    cur = conn.execute(
        """INSERT INTO ptm_site (uniprot_accession, position, residue, ptm_type_id,
                                  subtype, flanking_sequence, notes)
           VALUES (?, ?, ?, ?, ?, ?, ?)
           ON CONFLICT(uniprot_accession, position, ptm_type_id) DO UPDATE SET
             residue = COALESCE(excluded.residue, ptm_site.residue),
             subtype = COALESCE(excluded.subtype, ptm_site.subtype),
             flanking_sequence = COALESCE(excluded.flanking_sequence, ptm_site.flanking_sequence)
           RETURNING ptm_site_id""",
        (accession, position, residue, ptm_type_id, subtype, flanking_sequence, notes),
    )
    return cur.fetchone()[0]


def insert_evidence(conn, ptm_site_id, source_id, evidence_type,
                    evidence_code=None, method=None, pubmed_id=None,
                    confidence_score=None, raw_record=None):
    conn.execute(
        """INSERT INTO ptm_evidence (ptm_site_id, source_id, evidence_type, evidence_code,
                                      method, pubmed_id, confidence_score, date_retrieved, raw_record)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (ptm_site_id, source_id, evidence_type, evidence_code, method,
         pubmed_id, confidence_score, TODAY, json.dumps(raw_record) if raw_record else None),
    )


def upsert_protein_disease(conn, accession, doid, association_type=None,
                            source_id=None, pubmed_id=None, confidence=None, notes=None):
    conn.execute(
        """INSERT OR IGNORE INTO protein_disease
           (uniprot_accession, doid, association_type, source_id, pubmed_id, confidence, notes)
           VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (accession, doid, association_type, source_id, pubmed_id, confidence, notes),
    )


def upsert_ptm_site_disease(conn, ptm_site_id, doid, relationship=None,
                              direction=None, source_id=None, pubmed_id=None, notes=None):
    conn.execute(
        """INSERT OR IGNORE INTO ptm_site_disease
           (ptm_site_id, doid, relationship, direction, source_id, pubmed_id, notes)
           VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (ptm_site_id, doid, relationship, direction, source_id, pubmed_id, notes),
    )


def upsert_modifier(conn, modifier_id, modifier_type, name=None, mass=None,
                     formula=None, description=None):
    conn.execute(
        """INSERT OR IGNORE INTO modifier (modifier_id, modifier_type, name, mass, formula, description)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (modifier_id, modifier_type, name, mass, formula, description),
    )


def upsert_ptm_site_modifier(conn, ptm_site_id, modifier_id, source_id=None,
                               evidence_type=None, pubmed_id=None, notes=None):
    conn.execute(
        """INSERT OR IGNORE INTO ptm_site_modifier
           (ptm_site_id, modifier_id, source_id, evidence_type, pubmed_id, notes)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (ptm_site_id, modifier_id, source_id, evidence_type, pubmed_id, notes),
    )


def upsert_variant(conn, accession, position, ref_residue=None, alt_residue=None,
                    dbsnp_id=None, clinvar_id=None, clinical_significance=None, source_id=None):
    """Insert a variant; returns variant_id."""
    cur = conn.execute(
        """INSERT INTO variant (uniprot_accession, position, ref_residue, alt_residue,
                                 dbsnp_id, clinvar_id, clinical_significance, source_id)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?)
           RETURNING variant_id""",
        (accession, position, ref_residue, alt_residue, dbsnp_id, clinvar_id,
         clinical_significance, source_id),
    )
    return cur.fetchone()[0]


def upsert_tissue_expression(conn, gene_symbol, tissue, expression_tpm=None,
                              source_id=None, notes=None):
    conn.execute(
        """INSERT INTO tissue_expression (gene_symbol, tissue, expression_tpm, source_id, notes)
           VALUES (?, ?, ?, ?, ?)""",
        (gene_symbol, tissue, expression_tpm, source_id, notes),
    )
