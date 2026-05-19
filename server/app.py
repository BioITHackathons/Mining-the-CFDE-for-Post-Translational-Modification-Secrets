"""
PTM dashboard server.

Endpoints match the JSON files originally fetched by index.html. As the
DB-backed query layer (server/db.py) gains coverage, endpoints are
migrated from FileResponse-on-JSON to live SQLite queries one at a time.
"""

import os
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from server import db

REPO_ROOT = Path(__file__).resolve().parents[1]
DASHBOARD_DIR = REPO_ROOT / "dashboards" / "main" / "ptm-dashboard"

# Endpoints still backed by static JSON. Migrated routes are registered
# above this loop and removed from the dict as they move to db.py.
JSON_ENDPOINTS = {
    "sites": "sites.json",
    "quantitative": "quantitative.json",
    "enrichment": "enrichment.json",
    "site_evidence": "site_evidence.json",
    "glygen_data": "glygen_data.json",
}

app = FastAPI(title="CFDE PTM Dashboard API")

# gzip JSON payloads — /api/pqtl is ~45 MB raw / ~3 MB gzipped, /api/data
# similar ratio. Without compression, slow upload links + cloudflared
# combine to make the browser cancel mid-transfer.
app.add_middleware(GZipMiddleware, minimum_size=1000)

# Allow the dashboard to be hosted on a separate origin (a static bucket,
# Pages, localhost during dev) while the API runs on this server (typically
# exposed via a cloudflared tunnel). Override CORS_ORIGINS to restrict.
_cors_env = os.environ.get("CORS_ORIGINS", "*")
_cors_origins = [o.strip() for o in _cors_env.split(",") if o.strip()] or ["*"]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=False,
    allow_methods=["GET"],
    allow_headers=["*"],
)


# ── DB-backed routes ──────────────────────────────────────────────────

@app.get("/api/data")
async def api_data():
    return JSONResponse(db.fetch_data())


@app.get("/api/pqtl")
async def api_pqtl():
    return JSONResponse(db.fetch_pqtl())


# Sync def — blocking NCBI fetch runs in FastAPI's threadpool, so other
# requests can still be served while a slow protein-evidence call is in flight.
@app.get("/api/protein/{ac}/evidence")
def api_protein_evidence(ac: str):
    return JSONResponse(db.fetch_protein_evidence(ac))


@app.get("/api/protein/{ac}/quantitative")
def api_protein_quantitative(ac: str):
    return JSONResponse(db.fetch_protein_quantitative(ac))


@app.get("/api/tracks")
async def api_tracks():
    """Lightweight manifest only — full per-protein tracks load lazily."""
    return JSONResponse(db.fetch_tracks_manifest())


@app.get("/api/protein/{ac}/tracks")
def api_protein_tracks(ac: str):
    return JSONResponse(db.fetch_protein_tracks(ac))


@app.get("/api/pubmed/cache_stats")
async def api_pubmed_cache_stats():
    from server import pubmed
    return pubmed.cache_stats()


# ── Static-JSON routes (to be migrated) ───────────────────────────────

def _serve_json(filename: str) -> FileResponse:
    path = DASHBOARD_DIR / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"{filename} not found")
    return FileResponse(path, media_type="application/json")


for name, filename in JSON_ENDPOINTS.items():
    def _make_handler(fn: str):
        async def handler():
            return _serve_json(fn)
        return handler
    app.add_api_route(f"/api/{name}", _make_handler(filename), methods=["GET"])


@app.get("/api/health")
async def health():
    return {
        "status": "ok",
        "db_path": str(db.DB_PATH),
        "db_exists": db.DB_PATH.exists(),
        "live_endpoints": [
            "/api/data", "/api/pqtl", "/api/tracks",
            "/api/protein/{ac}/evidence",
            "/api/protein/{ac}/quantitative",
            "/api/protein/{ac}/tracks",
        ],
        "static_endpoints": [f"/api/{n}" for n in JSON_ENDPOINTS],
    }


app.mount("/", StaticFiles(directory=str(DASHBOARD_DIR), html=True), name="dashboard")
