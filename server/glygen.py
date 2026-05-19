"""
Cached fetcher for per-protein GlyGen data.

The pipeline loads GlyGen protein/site/glycan data into ptm_disease.db,
but it dropped per-variant disease associations on the floor — those came
back nested inside each `snv` object's `disease` list (BioMuta and friends).
We need that info to color the Track Viewer's variant diamonds by their
actual disease links, not the protein-level disease bag.

Rather than re-running the build, this module hits GlyGen's protein/detail
endpoint on demand and caches the response to disk (one JSON per protein).
First view of a protein takes ~500ms; subsequent loads are instant.
"""

from __future__ import annotations

import json
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = REPO_ROOT / "pipeline" / ".cache" / "glygen"


def _cache_path(uniprot_ac: str) -> Path:
    return CACHE_DIR / f"{uniprot_ac}.json"


def fetch_protein_detail(uniprot_ac: str, force: bool = False) -> dict:
    """Return GlyGen's full /protein/detail payload for `uniprot_ac`.

    Cached to disk forever — GlyGen entries are stable enough that we can
    treat them as immutable. Call with force=True to bust the cache.
    """
    path = _cache_path(uniprot_ac)
    if not force and path.exists():
        try:
            return json.loads(path.read_text())
        except Exception:
            pass  # corrupted cache; refetch

    url = f"https://api.glygen.org/protein/detail/{uniprot_ac}/"
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        data = r.json()
    except Exception:
        return {}

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return data


def snv_disease_map(uniprot_ac: str) -> dict[int, dict]:
    """Return {position: {diseases, variant_type, comment, sources, chr_id, chr_pos}}.

    Each variant in GlyGen carries its own (typically narrow) disease list
    sourced from BioMuta, COSMIC, etc. We key by start_pos because that's
    what the DB's `variant` table uses. Multiple alt residues at the same
    position get merged.
    """
    data = fetch_protein_detail(uniprot_ac)
    out: dict[int, dict] = {}
    for s in data.get("snv", []):
        pos = s.get("start_pos")
        if not isinstance(pos, int):
            continue
        diseases = []
        for d in s.get("disease", []):
            name = d.get("recommended_name", {}).get("name") or d.get("disease_id", "")
            if name and name not in diseases:
                diseases.append(name)
        sources = []
        chr_id = chr_pos = ""
        for ev in s.get("evidence", []):
            db = ev.get("database", "")
            if db and db not in sources:
                sources.append(db)
        # Some entries surface chr info in keywords or genomic_locus blocks
        for loc in s.get("chr_loc", []) or []:
            if not chr_id:
                chr_id = str(loc.get("chr_id", "") or "")
                chr_pos = str(loc.get("position", "") or "")
        comment = s.get("comment", "") or ""
        cl = comment.lower()
        vtype = "somatic" if "somatic" in cl else ("germline" if "germline" in cl else "")
        bucket = out.setdefault(pos, {
            "diseases": [], "variant_type": "", "comment": "",
            "sources": [], "chr_id": "", "chr_pos": "",
        })
        for d in diseases:
            if d not in bucket["diseases"]:
                bucket["diseases"].append(d)
        for src in sources:
            if src not in bucket["sources"]:
                bucket["sources"].append(src)
        if vtype and not bucket["variant_type"]:
            bucket["variant_type"] = vtype
        if comment and not bucket["comment"]:
            bucket["comment"] = comment
        if chr_id and not bucket["chr_id"]:
            bucket["chr_id"] = chr_id
            bucket["chr_pos"] = chr_pos
    return out
