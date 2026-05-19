-- ============================================================================
-- PTM-Disease Association Database Schema
-- Project: Mining Post-Translational Modifications in CFDE Resources
--          for Neurodegenerative Diseases
-- Version: 0.2
-- Engine:  SQLite
-- ============================================================================

PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

-- ============================================================================
-- 1. CORE ENTITIES
-- ============================================================================

-- Proteins (canonical human, keyed by UniProt accession)
CREATE TABLE protein (
    uniprot_accession   TEXT PRIMARY KEY,        -- e.g. "P05067"
    gene_symbol         TEXT NOT NULL,            -- HGNC symbol, e.g. "APP"
    protein_name        TEXT,                     -- full name
    organism            TEXT DEFAULT 'Homo sapiens',
    sequence_length     INTEGER,
    reviewed            INTEGER DEFAULT 1         -- 1 = Swiss-Prot, 0 = TrEMBL
);
CREATE INDEX idx_protein_gene ON protein(gene_symbol);

-- Diseases (keyed by Disease Ontology ID)
CREATE TABLE disease (
    doid                TEXT PRIMARY KEY,        -- e.g. "DOID:10652"
    name                TEXT NOT NULL,            -- "Alzheimer's disease"
    description         TEXT,
    category            TEXT                      -- "neurodegenerative", etc.
);

-- PTM types (keyed by GO term or PSI-MOD)
CREATE TABLE ptm_type (
    ptm_type_id         TEXT PRIMARY KEY,        -- GO or MOD accession
    name                TEXT NOT NULL,            -- "N-linked glycosylation"
    category            TEXT,                     -- "glycosylation", "phosphorylation", etc.
    description         TEXT
);

-- Data sources / CFDE DCCs
CREATE TABLE data_source (
    source_id           TEXT PRIMARY KEY,        -- e.g. "glygen", "idg", "lincs"
    name                TEXT NOT NULL,            -- "GlyGen"
    url                 TEXT,
    cfde_program        TEXT,                     -- CFDE program name if applicable
    description         TEXT
);

-- ============================================================================
-- 2. PTM SITES -- the central fact table
-- ============================================================================

-- One row = one (protein, position, ptm_type).
CREATE TABLE ptm_site (
    ptm_site_id         INTEGER PRIMARY KEY AUTOINCREMENT,
    uniprot_accession   TEXT NOT NULL REFERENCES protein(uniprot_accession),
    position            INTEGER NOT NULL,         -- 1-based residue position
    residue             TEXT,                      -- amino acid at that position
    ptm_type_id         TEXT NOT NULL REFERENCES ptm_type(ptm_type_id),
    subtype             TEXT,                      -- finer grain, e.g. "O-GlcNAcylation"
    flanking_sequence   TEXT,                      -- sequence window around the site
    notes               TEXT,

    UNIQUE(uniprot_accession, position, ptm_type_id)
);
CREATE INDEX idx_ptm_site_protein ON ptm_site(uniprot_accession);
CREATE INDEX idx_ptm_site_type    ON ptm_site(ptm_type_id);

-- ============================================================================
-- 3. EVIDENCE -- what supports each PTM site
-- ============================================================================

CREATE TABLE ptm_evidence (
    evidence_id         INTEGER PRIMARY KEY AUTOINCREMENT,
    ptm_site_id         INTEGER NOT NULL REFERENCES ptm_site(ptm_site_id),
    source_id           TEXT NOT NULL REFERENCES data_source(source_id),
    evidence_type       TEXT NOT NULL,            -- "experimental", "predicted", "text-mined", "inferred"
    evidence_code       TEXT,                     -- ECO term if available
    method              TEXT,                     -- e.g. "mass spectrometry", "motif prediction"
    pubmed_id           TEXT,
    confidence_score    REAL,                     -- 0-1 normalised
    date_retrieved      TEXT,                     -- ISO date
    raw_record          TEXT                      -- JSON blob of original source record
);
CREATE INDEX idx_evidence_site   ON ptm_evidence(ptm_site_id);
CREATE INDEX idx_evidence_source ON ptm_evidence(source_id);

-- ============================================================================
-- 4. DISEASE ASSOCIATIONS (two granularities)
-- ============================================================================

-- Protein-level disease links (many-to-many).
CREATE TABLE protein_disease (
    protein_disease_id  INTEGER PRIMARY KEY AUTOINCREMENT,
    uniprot_accession   TEXT NOT NULL REFERENCES protein(uniprot_accession),
    doid                TEXT NOT NULL REFERENCES disease(doid),
    association_type    TEXT,                     -- "causal", "genetic", "biomarker", "therapeutic_target"
    source_id           TEXT REFERENCES data_source(source_id),
    pubmed_id           TEXT,
    confidence          TEXT,                     -- "high", "medium", "low"
    notes               TEXT,

    UNIQUE(uniprot_accession, doid, association_type)
);
CREATE INDEX idx_pd_protein ON protein_disease(uniprot_accession);
CREATE INDEX idx_pd_disease ON protein_disease(doid);

-- Site-level disease links.
-- e.g. "hyperphosphorylation of tau at S396 in AD"
CREATE TABLE ptm_site_disease (
    ptm_site_disease_id INTEGER PRIMARY KEY AUTOINCREMENT,
    ptm_site_id         INTEGER NOT NULL REFERENCES ptm_site(ptm_site_id),
    doid                TEXT NOT NULL REFERENCES disease(doid),
    relationship        TEXT,                     -- "pathogenic", "protective", "biomarker", "altered_in"
    direction           TEXT,                     -- "increased", "decreased", "aberrant"
    source_id           TEXT REFERENCES data_source(source_id),
    pubmed_id           TEXT,
    notes               TEXT,

    UNIQUE(ptm_site_id, doid, relationship)
);
CREATE INDEX idx_psd_site    ON ptm_site_disease(ptm_site_id);
CREATE INDEX idx_psd_disease ON ptm_site_disease(doid);

-- ============================================================================
-- 5. MODIFIERS -- general entity attached at PTM sites
-- ============================================================================
-- Replaces glycan-specific tables. A "modifier" is any chemical group or
-- molecular entity attached to a protein at a PTM site:
--   glycan (GlyTouCan), phosphate group, ubiquitin, SUMO, acetyl, methyl, etc.
-- For simple modifications (phospho, acetyl) you may not need a modifier row
-- at all -- the ptm_type is sufficient. This table is for cases where the
-- attached entity has its own identity/structure (glycans, ubiquitin chains).

CREATE TABLE modifier (
    modifier_id         TEXT PRIMARY KEY,         -- canonical ID from source db
                                                  --   glycan: GlyTouCan acc ("G00055MO")
                                                  --   ubiquitin chain: "Ub-K48", "Ub-K63"
                                                  --   SUMO: "SUMO1", "SUMO2"
                                                  --   small molecule: ChEBI ID
    modifier_type       TEXT NOT NULL,            -- "glycan", "ubiquitin", "sumo", "small_molecule"
    name                TEXT,                     -- human-readable name
    mass                REAL,                     -- molecular mass (Da) if known
    formula             TEXT,                     -- molecular formula or composition string
    description         TEXT
);
CREATE INDEX idx_modifier_type ON modifier(modifier_type);

-- Junction: which modifiers are attached at which PTM sites.
CREATE TABLE ptm_site_modifier (
    ptm_site_modifier_id INTEGER PRIMARY KEY AUTOINCREMENT,
    ptm_site_id          INTEGER NOT NULL REFERENCES ptm_site(ptm_site_id),
    modifier_id          TEXT NOT NULL REFERENCES modifier(modifier_id),
    source_id            TEXT REFERENCES data_source(source_id),
    evidence_type        TEXT,                    -- "experimental", "predicted"
    pubmed_id            TEXT,
    notes                TEXT,

    UNIQUE(ptm_site_id, modifier_id)
);
CREATE INDEX idx_psm_site     ON ptm_site_modifier(ptm_site_id);
CREATE INDEX idx_psm_modifier ON ptm_site_modifier(modifier_id);

-- ============================================================================
-- 6. VARIANTS AT PTM SITES
-- ============================================================================

CREATE TABLE variant (
    variant_id          INTEGER PRIMARY KEY AUTOINCREMENT,
    uniprot_accession   TEXT NOT NULL REFERENCES protein(uniprot_accession),
    position            INTEGER NOT NULL,
    ref_residue         TEXT,
    alt_residue         TEXT,
    dbsnp_id            TEXT,                     -- e.g. "rs123456"
    clinvar_id          TEXT,
    clinical_significance TEXT,                   -- "pathogenic", "benign", "VUS"
    source_id           TEXT REFERENCES data_source(source_id)
);
CREATE INDEX idx_variant_protein ON variant(uniprot_accession, position);

CREATE TABLE variant_ptm_site (
    variant_ptm_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    variant_id          INTEGER NOT NULL REFERENCES variant(variant_id),
    ptm_site_id         INTEGER NOT NULL REFERENCES ptm_site(ptm_site_id),
    effect              TEXT,                     -- "abolishes_site", "creates_site", "adjacent", "unknown"
    distance            INTEGER DEFAULT 0,        -- residues from PTM site (0 = directly on it)
    notes               TEXT,

    UNIQUE(variant_id, ptm_site_id)
);

-- ============================================================================
-- 7. BIOMARKERS (source-agnostic)
-- ============================================================================

CREATE TABLE biomarker (
    biomarker_id        INTEGER PRIMARY KEY AUTOINCREMENT,
    name                TEXT NOT NULL,
    biomarker_type      TEXT,                     -- "protein", "modification", "metabolite"
    uniprot_accession   TEXT REFERENCES protein(uniprot_accession),
    ptm_site_id         INTEGER REFERENCES ptm_site(ptm_site_id),
    modifier_id         TEXT REFERENCES modifier(modifier_id),
    doid                TEXT REFERENCES disease(doid),
    specimen_type       TEXT,                     -- "serum", "CSF", "tissue"
    direction           TEXT,                     -- "elevated", "decreased"
    source_id           TEXT REFERENCES data_source(source_id),
    pubmed_id           TEXT
);
CREATE INDEX idx_biomarker_disease ON biomarker(doid);
CREATE INDEX idx_biomarker_protein ON biomarker(uniprot_accession);

-- ============================================================================
-- 8. CROSS-PTM CROSSTALK
-- ============================================================================

CREATE TABLE ptm_crosstalk (
    crosstalk_id        INTEGER PRIMARY KEY AUTOINCREMENT,
    ptm_site_id_a       INTEGER NOT NULL REFERENCES ptm_site(ptm_site_id),
    ptm_site_id_b       INTEGER NOT NULL REFERENCES ptm_site(ptm_site_id),
    relationship        TEXT,                     -- "competitive", "cooperative", "reciprocal", "proximity"
    distance_residues   INTEGER,
    source_id           TEXT REFERENCES data_source(source_id),
    pubmed_id           TEXT,
    notes               TEXT,

    CHECK(ptm_site_id_a < ptm_site_id_b)         -- avoid duplicate pairs
);

-- ============================================================================
-- 9. TISSUE EXPRESSION (GTEx / HuBMAP integration)
-- ============================================================================

CREATE TABLE tissue_expression (
    expression_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol         TEXT NOT NULL,
    tissue              TEXT NOT NULL,             -- UBERON term or GTEx label
    expression_tpm      REAL,
    source_id           TEXT REFERENCES data_source(source_id),
    notes               TEXT
);
CREATE INDEX idx_tissue_gene ON tissue_expression(gene_symbol);

-- ============================================================================
-- 10. VIEWS
-- ============================================================================

-- Flat summary: one row per (site, disease, evidence).
CREATE VIEW v_ptm_disease_summary AS
SELECT
    p.gene_symbol,
    p.uniprot_accession,
    p.protein_name,
    ps.position,
    ps.residue,
    pt.name            AS ptm_type,
    pt.category         AS ptm_category,
    ps.subtype,
    d.name              AS disease_name,
    d.doid,
    psd.relationship    AS disease_relationship,
    psd.direction       AS disease_direction,
    ds.name             AS data_source_name,
    e.evidence_type,
    e.confidence_score
FROM ptm_site ps
JOIN protein p         ON p.uniprot_accession = ps.uniprot_accession
JOIN ptm_type pt       ON pt.ptm_type_id = ps.ptm_type_id
LEFT JOIN ptm_site_disease psd ON psd.ptm_site_id = ps.ptm_site_id
LEFT JOIN disease d    ON d.doid = psd.doid
LEFT JOIN ptm_evidence e ON e.ptm_site_id = ps.ptm_site_id
LEFT JOIN data_source ds ON ds.source_id = e.source_id;

-- Crosstalk candidates: any two PTM types within 10 residues on the same protein.
CREATE VIEW v_crosstalk_candidates AS
SELECT
    p.gene_symbol,
    p.uniprot_accession,
    pt_a.name          AS ptm_type_a,
    a.position          AS position_a,
    a.residue           AS residue_a,
    pt_b.name          AS ptm_type_b,
    b.position          AS position_b,
    b.residue           AS residue_b,
    ABS(a.position - b.position) AS distance
FROM ptm_site a
JOIN ptm_site b    ON a.uniprot_accession = b.uniprot_accession
                   AND a.ptm_site_id < b.ptm_site_id
                   AND a.ptm_type_id != b.ptm_type_id
JOIN protein p     ON p.uniprot_accession = a.uniprot_accession
JOIN ptm_type pt_a ON pt_a.ptm_type_id = a.ptm_type_id
JOIN ptm_type pt_b ON pt_b.ptm_type_id = b.ptm_type_id
WHERE ABS(a.position - b.position) <= 10
ORDER BY p.gene_symbol, distance;

-- Sites with attached modifiers (e.g. specific glycan structures, ubiquitin chains).
CREATE VIEW v_site_modifiers AS
SELECT
    p.gene_symbol,
    p.uniprot_accession,
    ps.position,
    ps.residue,
    pt.name            AS ptm_type,
    m.modifier_id,
    m.modifier_type,
    m.name             AS modifier_name,
    m.mass,
    ds.name            AS data_source_name
FROM ptm_site_modifier psm
JOIN ptm_site ps   ON ps.ptm_site_id = psm.ptm_site_id
JOIN protein p     ON p.uniprot_accession = ps.uniprot_accession
JOIN ptm_type pt   ON pt.ptm_type_id = ps.ptm_type_id
JOIN modifier m    ON m.modifier_id = psm.modifier_id
LEFT JOIN data_source ds ON ds.source_id = psm.source_id;
