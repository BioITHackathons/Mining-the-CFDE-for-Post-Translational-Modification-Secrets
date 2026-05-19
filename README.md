# Mining the CFDE for Post-Translational Modification Secrets

A live, queryable PTM-disease knowledge base assembled from GlyGen, MoTrPAC,
and CFDE Knowledge Centers, surfaced through a FastAPI server with real-time
PubMed mining and a browser dashboard.

## Repository layout

```
pipeline/                         Build the SQLite DB from upstream sources
  build_data.py                   Resumable, staged loader (seed → glygen → motrpac → kc)
  ptm_disease_db_schema.sql       Schema
  loaders/                        Per-source loader helpers
  ptm_disease.db                  Built artifact (gitignored)
  .cache/                         Runtime caches (PubMed abstracts; gitignored)
  .checkpoints/                   Stage checkpoints (gitignored)

server/                           FastAPI app + DB query layer
  app.py                          Routes (+ CORS)
  db.py                           Per-protein assembly: data, pqtl, tracks, evidence, quant
  pubmed.py                       NCBI E-utils fetch + SQLite cache + regex NLP

dashboards/
  main/
    motrpac-dashboard/            MoTrPAC PTM Exercise Response Explorer
    ptm-dashboard/                CFDE/GlyGen PTM Disease-Tissue Explorer
      index.html                  Pure HTML/JS — talks to /api/* (configurable origin)
      legacy/                     Original offline build scripts (superseded by server/)
  other/                          Earlier prototype dashboards

serve.sh                          Local launcher (FastAPI + datasette + cloudflared)
pyproject.toml                    uv project
```

## Quickstart

```bash
uv sync                                  # install deps
uv run python pipeline/build_data.py     # build the SQLite DB (long; resumable)
./serve.sh                               # FastAPI + datasette + cloudflared tunnel
```

`./serve.sh` starts three processes:

- **FastAPI** on `http://127.0.0.1:8000` — serves the dashboard at `/` and the API under `/api/*`
- **datasette** on `http://127.0.0.1:8001` — ad-hoc SQL UI over the same DB
- **cloudflared** — opens a public tunnel and prints a `https://<random>.trycloudflare.com` URL pointing at the FastAPI port

Override via env: `APP_PORT=9000 DB=path/to.db ./serve.sh`.

## API endpoints

All endpoints return JSON. The dashboard's startup load uses the first six.

| Path | Backend | Notes |
|---|---|---|
| `GET /api/data` | live DB | `{summary, proteins[]}` — every protein in the DB; NLP fields merged from curated `glygen_data.json` where available |
| `GET /api/pqtl` | live DB | `{summary, overlaps[], protein_summaries{}}` — variant ↔ PTM-site joins with disease links |
| `GET /api/tracks` | live DB | Lightweight manifest: protein names + site counts; per-protein detail loads lazily |
| `GET /api/protein/{ac}/tracks` | live DB + cached NCBI | Full single-protein track payload (sites, conditions, scores, snvs) |
| `GET /api/protein/{ac}/evidence` | live DB + cached NCBI | Per-site PubMed-mined conditions/tissues/methods |
| `GET /api/protein/{ac}/quantitative` | live DB + cached NCBI | Sentence-level direction/fold-change hits |
| `GET /api/sites` `quantitative` `enrichment` | static JSON | Curated 30-protein snapshots (still used as defaults) |
| `GET /api/pubmed/cache_stats` | meta | PubMed abstract cache hits |
| `GET /api/health` | meta | DB path, live vs static endpoints |

PubMed abstracts are fetched once from NCBI E-utils and cached in
`pipeline/.cache/pubmed_cache.db` (SQLite, WAL mode). Cold per-protein
fetches take ~200–700ms; warm refetches are ~150ms.

## Deployment: host the dashboard separately from the API

The dashboard is a single static HTML file. You can serve it from a GCS
bucket, S3, GitHub Pages, etc., and point it at an API running anywhere —
typically your laptop, exposed through cloudflared.

1. Start the API locally with `./serve.sh`. Copy the `trycloudflare.com`
   URL it prints.
2. Tell the dashboard where to find the API. Two options:

   **Edit the file before uploading:**
   ```html
   <script>
     window.API_BASE = 'https://random-name.trycloudflare.com';
   </script>
   ```
   (See the top of `dashboards/main/ptm-dashboard/index.html`.)

   **Override at runtime via URL param** (handy if the tunnel rotates):
   ```
   https://your-bucket.example/index.html?api=https://random-name.trycloudflare.com
   ```

3. Upload `index.html` (and any sibling JSONs it still references) to your
   bucket. No other build step.

CORS is open by default (`*`). Restrict to your bucket origin by setting
`CORS_ORIGINS=https://your-bucket.example` before launching the server.

## Database build

The build is split into independently resumable stages. State is tracked
in `pipeline/.checkpoints/` and in the DB itself.

```bash
uv run python pipeline/build_data.py                  # all stages
uv run python pipeline/build_data.py --stage seed     # ref proteins/diseases
uv run python pipeline/build_data.py --stage glygen   # all human glycoproteins (long)
uv run python pipeline/build_data.py --stage motrpac  # phospho/acetyl/ubiquityl
uv run python pipeline/build_data.py --stage kc       # GTEx / Kids First / IDG
uv run python pipeline/build_data.py --stage summary  # print stats
uv run python pipeline/build_data.py --reset          # delete DB + checkpoints
```

Re-running any stage skips items it's already processed.

## Architecture notes

- The legacy offline build scripts in `dashboards/main/ptm-dashboard/legacy/`
  produced the JSON files the dashboard originally loaded. They've been
  replaced by live `/api/*` endpoints — see that directory's `README.md`
  for the full mapping.
- The dashboard never queries the DB directly — it goes through `/api/*`.
- The PubMed mining (regex term-dictionaries + sentence-level direction
  extraction) lives in `server/pubmed.py`. Extending the vocabularies
  there affects every live request.
- `glygen_data.json` is still referenced as the source of NLP-derived
  fields (`expression_disease`, `ptm_annotations_*`, `publication_count`)
  for the curated 30 proteins. Other proteins get those fields empty
  until a `stage_pubmed_curated` pipeline stage is added.

## License

See [LICENSE](LICENSE).
