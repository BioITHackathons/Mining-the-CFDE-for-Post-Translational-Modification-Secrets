# Legacy build scripts

These scripts produced the JSON files the dashboard originally consumed.
They've been superseded by live API endpoints served by `server/app.py`
(see repo root). The scripts remain here for reference — to inspect the
original NLP rules and scoring formulas, or to rebuild the static JSONs
if you need an offline snapshot.

| Script | Output JSON | Replaced by |
|---|---|---|
| `fetch_glygen.py` | `glygen_data.json` | `pipeline/build_data.py` (GlyGen stage writes to `ptm_disease.db`) + `/api/data` |
| `fetch_site_evidence.py` | `site_evidence.json` | `/api/protein/{ac}/evidence` — real-time PubMed mining with SQLite cache |
| `extract_quantitative.py` | `quantitative.json` | `/api/protein/{ac}/quantitative` |
| `build_tracks.py` | `tracks.json` | `/api/tracks` (manifest) + `/api/protein/{ac}/tracks` |
| `build_scores.py` | (in-place merge into `tracks.json`) | Inline in `fetch_protein_tracks` — same scoring formula, including the protein-level fallback |
| `build_enrichment.py` | `enrichment.json` | Still served from the static JSON. Port to a pipeline stage when ready. |
| `build_pqtl.py` | `pqtl.json` | `/api/pqtl` — DB-side JOIN on `variant_ptm_site` + `protein_disease` |

## Hardcoded inputs that don't apply anymore

`fetch_glygen.py:11-42` and `fetch_site_evidence.py:288-295` both hardcoded
30 UniProt accessions for the curated neurodegeneration set. The live API
draws from every protein in `ptm_disease.db` (currently 2400+) — the
hardcoded list was a hackathon-day expedient, not a design constraint.

## NLP rules (lifted from these scripts)

The regex term dictionaries and `mine_sentence` / `mine_text` helpers
were lifted verbatim into `server/pubmed.py`. If you want to extend the
vocabularies (more diseases, more tissues, more direction words), edit
those server-side dicts — they affect every live-mining request.

## Re-generating static JSONs

If you really need to rebuild the static files (e.g., to bundle them with
the dashboard for a snapshot release), the original order was:

```
python fetch_glygen.py
python fetch_site_evidence.py
python extract_quantitative.py
python build_tracks.py
python build_scores.py
python build_enrichment.py
python build_pqtl.py
```

Each script reads/writes hardcoded `/app/workspace/...` paths that won't
exist locally — you'll need to patch the file paths to your repo root
before running them.
