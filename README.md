# Mining the CFDE for Post-Translational Modification Secrets

## Repository Structure

```
dashboards/
  main/
    motrpac-dashboard/   # MoTrPAC PTM Exercise Response Explorer
    ptm-dashboard/       # CFDE/GlyGen PTM Disease-Tissue Explorer
  other/
    dashboard/           # Original network/pathway dashboard
    ptm-3d-viewer/       # Protein 3D Structure Viewer (AlphaFold)
    ptm-cooccurrence/    # PTM Co-occurrence Explorer
    ptm-gaps-explorer/   # PTM Research Gaps Explorer
    ptm-strategy/        # PTM Strategy Dashboard
    ptm-pathway-network/ # PTM Pathway Network Viewer

pipeline/
  loaders/               # Database loader scripts
  build_data.py          # Main build script
  ptm_disease.db         # Pre-built SQLite database
  ptm_disease_db_schema.sql  # Database schema
```

## Dashboards

### Main

- **MoTrPAC Dashboard** (`dashboards/main/motrpac-dashboard/`) — Interactive explorer for exercise-induced PTM changes from MoTrPAC data.
- **PTM Dashboard** (`dashboards/main/ptm-dashboard/`) — Unified CFDE/GlyGen PTM disease-tissue explorer with enrichment analysis, pQTL data, and site-level evidence.

### Other

- **Dashboard** — Original pathway network dashboard.
- **3D Viewer** — AlphaFold-based protein structure viewer with PTM site annotations.
- **Co-occurrence** — PTM co-occurrence and crosstalk analysis.
- **Gaps Explorer** — Research gaps analysis across PTM types and tissues.
- **Strategy** — Strategic overview of PTM research priorities.
- **Pathway Network** — PTM pathway network visualization.

## Pipeline

The `pipeline/` directory contains scripts to build the PTM disease database from CFDE, GlyGen, and MoTrPAC sources.

```bash
# Build the database
cd pipeline
python build_data.py
```

## License

See [LICENSE](LICENSE).
