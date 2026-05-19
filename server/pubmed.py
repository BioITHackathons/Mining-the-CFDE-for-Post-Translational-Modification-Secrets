"""
PubMed abstract fetch + lightweight regex NLP.

Lifts the term dictionaries and `mine_text` regex logic from
dashboards/main/ptm-dashboard/fetch_site_evidence.py so the live API can
do the same extraction at request time. Adds a small SQLite cache so a
given PMID is only fetched from NCBI once across all sessions.
"""

from __future__ import annotations

import re
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_BATCH = 100
NCBI_SLEEP = 0.4

REPO_ROOT = Path(__file__).resolve().parents[1]
CACHE_PATH = REPO_ROOT / "pipeline" / ".cache" / "pubmed_cache.db"


# ── Term dictionaries (lifted from fetch_site_evidence.py) ─────────────

DISEASE_TERMS = {
    r"\balzheimer": "Alzheimer's disease",
    r"\bparkinson": "Parkinson's disease",
    r"\bhuntington": "Huntington's disease",
    r"\bamyotrophic lateral sclerosis\b": "ALS",
    r"\b[Aa][Ll][Ss]\b": "ALS",
    r"\bprion\b": "Prion disease",
    r"\bcreutzfeldt": "Creutzfeldt-Jakob disease",
    r"\bdementia\b": "Dementia",
    r"\btauopath": "Tauopathy",
    r"\bsynucleinopath": "Synucleinopathy",
    r"\blewy bod": "Lewy body disease",
    r"\bfrontotemporal": "Frontotemporal dementia",
    r"\bmultiple sclerosis": "Multiple sclerosis",
    r"\bepilep": "Epilepsy",
    r"\bstroke\b": "Stroke",
    r"\bischemi": "Ischemia",
    r"\bneurodegenerat": "Neurodegeneration",
    r"\bneuroinflam": "Neuroinflammation",
    r"\bglioblastoma": "Glioblastoma",
    r"\bglioma\b": "Glioma",
    r"\bmelanoma\b": "Melanoma",
    r"\bbreast cancer": "Breast cancer",
    r"\blung cancer": "Lung cancer",
    r"\bcolorectal cancer": "Colorectal cancer",
    r"\bpancreatic cancer": "Pancreatic cancer",
    r"\bprostate cancer": "Prostate cancer",
    r"\bhepatocellular carcinoma": "Hepatocellular carcinoma",
    r"\bovarian cancer": "Ovarian cancer",
    r"\bgastric cancer": "Gastric cancer",
    r"\bcancer\b": "Cancer (general)",
    r"\btumor\b": "Tumor",
    r"\bcarcinoma\b": "Carcinoma",
    r"\bdiabet": "Diabetes",
    r"\binsulin resist": "Insulin resistance",
    r"\bobes": "Obesity",
    r"\batheros": "Atherosclerosis",
    r"\bcardiovasc": "Cardiovascular disease",
    r"\bheart fail": "Heart failure",
    r"\bfibrosis\b": "Fibrosis",
    r"\binflam": "Inflammation",
    r"\bautoimmun": "Autoimmune disease",
    r"\brheumatoid": "Rheumatoid arthritis",
    r"\blupus\b": "Lupus",
    r"\bsepsis\b": "Sepsis",
    r"\binfect": "Infection",
    r"\bviral\b": "Viral infection",
    r"\bCOVID": "COVID-19",
    r"\bSARS": "SARS-CoV-2",
    r"\bamyloid": "Amyloidosis",
    r"\bneurotox": "Neurotoxicity",
    r"\bapoptosis\b": "Apoptosis",
    r"\bdown syndrome": "Down syndrome",
    r"\bniemann.pick": "Niemann-Pick disease",
    r"\bnormal\b": "Normal/healthy",
    r"\bhealthy\b": "Normal/healthy",
    r"\bcontrol\b": "Control condition",
    r"\baging\b": "Aging",
    r"\baged\b": "Aging",
}

TISSUE_TERMS = {
    r"\bbrain\b": "Brain",
    r"\bcerebr": "Cerebral cortex",
    r"\bcortex\b": "Cortex",
    r"\bcortical\b": "Cortex",
    r"\bhippocamp": "Hippocampus",
    r"\bsubstantia nigra": "Substantia nigra",
    r"\bstriatum\b": "Striatum",
    r"\bthalamus\b": "Thalamus",
    r"\bhypothalam": "Hypothalamus",
    r"\bcerebellum\b": "Cerebellum",
    r"\bamygdala\b": "Amygdala",
    r"\btemporal lobe": "Temporal lobe",
    r"\bfrontal cortex": "Frontal cortex",
    r"\bprefrontal": "Prefrontal cortex",
    r"\bdorsolateral": "Dorsolateral prefrontal cortex",
    r"\bspinal cord": "Spinal cord",
    r"\bneuron": "Neuron",
    r"\bastrocyte": "Astrocyte",
    r"\bmicroglia": "Microglia",
    r"\boligodendroc": "Oligodendrocyte",
    r"\bglia\b": "Glia",
    r"\bsynap": "Synapse",
    r"\baxon": "Axon",
    r"\bdendrit": "Dendrite",
    r"\bdopaminerg": "Dopaminergic neuron",
    r"\bcholinerg": "Cholinergic neuron",
    r"\bGABAerg": "GABAergic neuron",
    r"\bglutamaterg": "Glutamatergic neuron",
    r"\bliver\b": "Liver",
    r"\bhepat": "Liver/Hepatocyte",
    r"\bkidney\b": "Kidney",
    r"\brenal\b": "Kidney",
    r"\bheart\b": "Heart",
    r"\bcardiac\b": "Heart",
    r"\blung\b": "Lung",
    r"\bpulmonar": "Lung",
    r"\bpancrea": "Pancreas",
    r"\bintestin": "Intestine",
    r"\bcolon\b": "Colon",
    r"\bstomach\b": "Stomach",
    r"\bprostate\b": "Prostate",
    r"\bbreast\b": "Breast",
    r"\bovary\b": "Ovary",
    r"\btestes\b": "Testes",
    r"\bplasma\b": "Blood plasma",
    r"\bserum\b": "Blood serum",
    r"\bblood\b": "Blood",
    r"\berythrocyte": "Erythrocyte",
    r"\bplatelet": "Platelet",
    r"\bCSF\b": "Cerebrospinal fluid",
    r"\bcerebrospinal fluid": "Cerebrospinal fluid",
    r"\bendotheli": "Endothelium",
    r"\bepitheli": "Epithelium",
    r"\bfibroblast": "Fibroblast",
    r"\bmacrophage": "Macrophage",
    r"\bmonocyte": "Monocyte",
    r"\bT cell": "T cell",
    r"\bB cell": "B cell",
    r"\bNK cell": "NK cell",
    r"\bskeletal muscle": "Skeletal muscle",
    r"\bmuscle\b": "Muscle",
    r"\badipose": "Adipose tissue",
    r"\bbone\b": "Bone",
    r"\bskin\b": "Skin",
    r"\bpost.mortem": "Post-mortem tissue",
    r"\bautops": "Autopsy tissue",
    r"\bbiopsy\b": "Biopsy",
    r"\bcell line\b": "Cell line",
    r"\bHEK.?293": "HEK293 cells",
    r"\bHeLa\b": "HeLa cells",
    r"\bSH.SY5Y": "SH-SY5Y cells",
    r"\bin vitro\b": "In vitro",
    r"\bin vivo\b": "In vivo",
    r"\bmouse\b": "Mouse",
    r"\brat\b": "Rat",
    r"\bhuman\b": "Human",
}

EXPERIMENTAL_TERMS = {
    r"\bmass spectrom": "Mass spectrometry",
    r"\bLC.MS": "LC-MS/MS",
    r"\bproteom": "Proteomics",
    r"\bglycoprot": "Glycoproteomics",
    r"\bphosphoproteo": "Phosphoproteomics",
    r"\bimmunoprecip": "Immunoprecipitation",
    r"\bwestern blot": "Western blot",
    r"\bELISA\b": "ELISA",
    r"\blectin\b": "Lectin-based assay",
    r"\bcrystallog": "X-ray crystallography",
    r"\bcryo.EM": "Cryo-EM",
    r"\bNMR\b": "NMR",
    r"\bsite.directed mutag": "Site-directed mutagenesis",
    r"\bknockout\b": "Knockout",
    r"\bknockdown\b": "Knockdown",
    r"\btransgenic\b": "Transgenic model",
    r"\bTMT\b": "TMT labeling",
    r"\biTRAQ\b": "iTRAQ labeling",
    r"\blabel.free": "Label-free quantification",
    r"\benrichment\b": "Enrichment",
    r"\bHILIC\b": "HILIC enrichment",
    r"\bTiO2\b": "TiO2 enrichment",
    r"\bIMAC\b": "IMAC enrichment",
}

_COMPILED = {
    "disease": [(re.compile(p, re.IGNORECASE), n) for p, n in DISEASE_TERMS.items()],
    "tissue": [(re.compile(p, re.IGNORECASE), n) for p, n in TISSUE_TERMS.items()],
    "method": [(re.compile(p, re.IGNORECASE), n) for p, n in EXPERIMENTAL_TERMS.items()],
}


# ── Quantitative mining (lifted from extract_quantitative.py) ─────────

SITE_PAT = re.compile(
    r"(?:"
    r"(?:Ser|Thr|Tyr|Lys|Asn|Cys|Arg|His|Gly|Ala|Val|Leu|Ile|Pro|Phe|Met|Trp|Asp|Glu|Gln)"
    r"[\s\-]?(\d{1,4})"
    r"|"
    r"(?:p?[STYNKCRHDEGAQWFLMIPV])(\d{2,4})\b"
    r"|"
    r"(?:position|residue|site)\s+(\d{2,4})"
    r")",
    re.IGNORECASE,
)

FOLD_PAT = re.compile(r"(\d+\.?\d*)\s*-?\s*fold", re.IGNORECASE)
PCT_PAT = re.compile(r"(\d+\.?\d*)\s*%", re.IGNORECASE)
RATIO_PAT = re.compile(r"ratio\s+(?:of\s+)?(\d+\.?\d*)", re.IGNORECASE)

UP_WORDS = (
    "increased", "upregulated", "up-regulated", "elevated", "enhanced",
    "higher", "more prevalent", "accumulation", "hyperphosphorylat",
    "overexpress", "gain", "promotes", "activates", "stimulates",
    "induces", "augment", "amplif", "enrich", "abundant",
)
DOWN_WORDS = (
    "decreased", "downregulated", "down-regulated", "reduced", "lower",
    "loss", "inhibits", "blocks", "prevents", "impaired", "deficient",
    "absent", "abolishes", "disrupts", "attenuate", "diminish",
    "deplet", "suppress",
)

# Narrower per-sentence vocabularies — kept distinct from the broader
# DISEASE_TERMS/TISSUE_TERMS used by mine_text(), so the quantitative
# pass stays faithful to extract_quantitative.py.
DISEASE_MAP = {
    r"\balzheimer": "Alzheimer's disease",
    r"\bparkinson": "Parkinson's disease",
    r"\bhuntington": "Huntington's disease",
    r"\bALS\b": "ALS",
    r"\bamyotrophic lateral sclerosis": "ALS",
    r"\bprion\b": "Prion disease",
    r"\bdementia\b": "Dementia",
    r"\btauopath": "Tauopathy",
    r"\bsynucleinopath": "Synucleinopathy",
    r"\blewy bod": "Lewy body disease",
    r"\bfrontotemporal": "Frontotemporal dementia",
    r"\bmultiple sclerosis": "Multiple sclerosis",
    r"\bepilep": "Epilepsy",
    r"\bstroke\b": "Stroke",
    r"\bischemi": "Ischemia",
    r"\bneurodegenerat": "Neurodegeneration",
    r"\bneuroinflam": "Neuroinflammation",
    r"\bglioblastoma": "Glioblastoma",
    r"\bglioma\b": "Glioma",
    r"\bmelanoma\b": "Melanoma",
    r"\bbreast cancer": "Breast cancer",
    r"\blung cancer": "Lung cancer",
    r"\bcolorectal": "Colorectal cancer",
    r"\bpancreatic cancer": "Pancreatic cancer",
    r"\bprostate cancer": "Prostate cancer",
    r"\bhepatocellular": "Hepatocellular carcinoma",
    r"\bovarian cancer": "Ovarian cancer",
    r"\bgastric cancer": "Gastric cancer",
    r"\bcancer\b": "Cancer",
    r"\btumor\b": "Tumor",
    r"\bcarcinoma\b": "Carcinoma",
    r"\bdiabet": "Diabetes",
    r"\bobes": "Obesity",
    r"\batheros": "Atherosclerosis",
    r"\bcardiovasc": "Cardiovascular disease",
    r"\bfibrosis\b": "Fibrosis",
    r"\binflam": "Inflammation",
    r"\bautoimmun": "Autoimmune disease",
    r"\brheumatoid": "Rheumatoid arthritis",
    r"\bsepsis\b": "Sepsis",
    r"\bamyloid": "Amyloidosis",
    r"\baging\b": "Aging",
    r"\bnormal\b": "Normal/healthy",
    r"\bhealthy\b": "Normal/healthy",
    r"\bcontrol\b": "Control",
    r"\bdown syndrome": "Down syndrome",
    r"\bniemann.pick": "Niemann-Pick disease",
    r"\bsclero": "Sclerosis",
    r"\bcerebral amyloid": "Cerebral amyloid angiopathy",
    r"\bhypoxia": "Hypoxia",
    r"\boxidative stress": "Oxidative stress",
    r"\bER stress": "ER stress",
    r"\bapoptosis\b": "Apoptosis",
}

TISSUE_MAP = {
    r"\bbrain\b": "Brain",
    r"\bhippocamp": "Hippocampus",
    r"\bcortex\b": "Cortex",
    r"\bsubstantia nigra": "Substantia nigra",
    r"\bstriatum": "Striatum",
    r"\bcerebellum": "Cerebellum",
    r"\bneuron": "Neuron",
    r"\bastrocyte": "Astrocyte",
    r"\bmicroglia": "Microglia",
    r"\bsynap": "Synapse",
    r"\bdopaminerg": "Dopaminergic neuron",
    r"\bliver\b": "Liver",
    r"\bkidney\b": "Kidney",
    r"\bheart\b": "Heart",
    r"\blung\b": "Lung",
    r"\bpancrea": "Pancreas",
    r"\bplasma\b": "Plasma",
    r"\bserum\b": "Serum",
    r"\bblood\b": "Blood",
    r"\bCSF\b": "CSF",
    r"\bcerebrospinal": "CSF",
    r"\bmuscle\b": "Muscle",
    r"\bpost.mortem": "Post-mortem tissue",
    r"\bcell line": "Cell line",
    r"\bhuman\b": "Human",
    r"\bmouse\b": "Mouse",
    r"\brat\b": "Rat",
    r"\bin vitro": "In vitro",
    r"\bin vivo": "In vivo",
}

PTM_MAP = {
    r"\bphosphorylat": "Phosphorylation",
    r"\bglycosylat": "Glycosylation",
    r"\bn-glyco": "N-glycosylation",
    r"\bo-glyco": "O-glycosylation",
    r"\bubiquitin": "Ubiquitination",
    r"\bsumoylat": "SUMOylation",
    r"\bacetylat": "Acetylation",
    r"\bmethylat": "Methylation",
    r"\bpalmitoyl": "Palmitoylation",
    r"\bsulfat": "Sulfation",
    r"\bnitrosylat": "S-nitrosylation",
    r"\bsialylat": "Sialylation",
    r"\bfucosylat": "Fucosylation",
    r"\bglycation": "Glycation",
    r"\bcleavage": "Proteolytic cleavage",
}

_DISEASE_MAP_C = [(re.compile(p, re.IGNORECASE), n) for p, n in DISEASE_MAP.items()]
_TISSUE_MAP_C = [(re.compile(p, re.IGNORECASE), n) for p, n in TISSUE_MAP.items()]
_PTM_MAP_C = [(re.compile(p, re.IGNORECASE), n) for p, n in PTM_MAP.items()]


def extract_direction(text: str) -> str | None:
    tl = text.lower()
    for w in UP_WORDS:
        if w in tl:
            return "up"
    for w in DOWN_WORDS:
        if w in tl:
            return "down"
    return None


def extract_quant(text: str) -> dict | None:
    """Pick the first fold-change / percent / ratio mention in the sentence."""
    fm = FOLD_PAT.search(text)
    if fm:
        return {"type": "fold_change", "value": float(fm.group(1)), "raw": fm.group(0)}
    pm = PCT_PAT.search(text)
    if pm:
        return {"type": "percent", "value": float(pm.group(1)), "raw": pm.group(0)}
    rm = RATIO_PAT.search(text)
    if rm:
        return {"type": "ratio", "value": float(rm.group(1)), "raw": rm.group(0)}
    return None


def _match_compiled(text: str, compiled) -> list[str]:
    out = set()
    for pat, name in compiled:
        if pat.search(text):
            out.add(name)
    return sorted(out)


def mine_sentence(sent: str, protein_ac: str, protein_name: str, gene_name: str) -> dict | None:
    """Look for a site-level quantitative claim in one sentence.

    Returns a dict if the sentence has a direction word + at least one of
    {site, disease, tissue, ptm} mention; otherwise None.
    """
    sites: list[str] = []
    for m in SITE_PAT.finditer(sent):
        pos = m.group(1) or m.group(2) or m.group(3)
        if not pos:
            continue
        pre = sent[max(0, m.start() - 4):m.start()]
        aa_m = re.search(r"(Ser|Thr|Tyr|Lys|Asn|Cys|Arg|Glu|Asp|Gln|His|[STYNKC])\s*$", pre, re.IGNORECASE)
        aa = aa_m.group(1) if aa_m else ""
        sites.append(f"{aa}{pos}".strip())

    diseases = _match_compiled(sent, _DISEASE_MAP_C)
    tissues = _match_compiled(sent, _TISSUE_MAP_C)
    ptms = _match_compiled(sent, _PTM_MAP_C)
    direction = extract_direction(sent)
    quant = extract_quant(sent)

    if direction and (diseases or tissues or ptms or sites):
        return {
            "protein_ac": protein_ac,
            "protein_name": protein_name,
            "gene_name": gene_name,
            "sites": sites,
            "diseases": diseases,
            "tissues": tissues,
            "ptm_types": ptms,
            "direction": direction,
            "quantitative": quant,
            "sentence": sent.strip(),
        }
    return None


def split_sentences(text: str) -> list[str]:
    """Sentence split matching extract_quantitative.py (look-behind on . ! ?)."""
    if not text:
        return []
    return re.split(r"(?<=[.!?])\s+", text)


def mine_text(text: str) -> tuple[list[str], list[str], list[str]]:
    """Return (diseases, tissues, methods) found in text."""
    if not text:
        return [], [], []
    diseases, tissues, methods = set(), set(), set()
    for pat, name in _COMPILED["disease"]:
        if pat.search(text):
            diseases.add(name)
    for pat, name in _COMPILED["tissue"]:
        if pat.search(text):
            tissues.add(name)
    for pat, name in _COMPILED["method"]:
        if pat.search(text):
            methods.add(name)
    return sorted(diseases), sorted(tissues), sorted(methods)


# ── SQLite abstract cache ─────────────────────────────────────────────

_CACHE_INIT_DONE = False


def _cache_conn() -> sqlite3.Connection:
    global _CACHE_INIT_DONE
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(CACHE_PATH))
    conn.execute("PRAGMA journal_mode = WAL")
    if not _CACHE_INIT_DONE:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS pubmed_abstract (
                pmid TEXT PRIMARY KEY,
                title TEXT,
                abstract TEXT,
                fetched_at TEXT
            )
        """)
        _CACHE_INIT_DONE = True
    return conn


def _cache_get(pmids: list[str]) -> dict[str, dict]:
    if not pmids:
        return {}
    conn = _cache_conn()
    placeholders = ",".join("?" * len(pmids))
    rows = conn.execute(
        f"SELECT pmid, title, abstract FROM pubmed_abstract WHERE pmid IN ({placeholders})",
        pmids,
    ).fetchall()
    conn.close()
    return {r[0]: {"title": r[1] or "", "abstract": r[2] or ""} for r in rows}


def _cache_put(records: dict[str, dict]) -> None:
    if not records:
        return
    now = datetime.now(timezone.utc).isoformat()
    conn = _cache_conn()
    conn.executemany(
        "INSERT OR REPLACE INTO pubmed_abstract (pmid, title, abstract, fetched_at) VALUES (?, ?, ?, ?)",
        [(pmid, rec.get("title", ""), rec.get("abstract", ""), now) for pmid, rec in records.items()],
    )
    conn.commit()
    conn.close()


# ── NCBI E-utils fetch ────────────────────────────────────────────────

def _parse_pubmed_xml(xml_text: str, requested: list[str]) -> dict[str, dict]:
    """Extract per-PMID {title, abstract} from a PubmedArticleSet XML blob.

    PubmedArticle blocks aren't always in the same order as the request,
    and multi-PMID errors return partial XML. We scan article-by-article.
    """
    out: dict[str, dict] = {}
    articles = re.findall(r"<PubmedArticle>(.*?)</PubmedArticle>", xml_text, re.DOTALL)
    for art in articles:
        m_pmid = re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
        if not m_pmid:
            continue
        pmid = m_pmid.group(1)
        # ArticleTitle is required
        m_title = re.search(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", art, re.DOTALL)
        title = re.sub(r"<[^>]+>", "", m_title.group(1)).strip() if m_title else ""
        # Abstract can have multiple <AbstractText> blocks (structured abstracts)
        abs_parts = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, re.DOTALL)
        abstract = " ".join(re.sub(r"<[^>]+>", "", p).strip() for p in abs_parts).strip()
        out[pmid] = {"title": title, "abstract": abstract}
    # Mark any requested PMIDs not found in the response so we don't refetch them
    for pmid in requested:
        out.setdefault(pmid, {"title": "", "abstract": ""})
    return out


def fetch_abstracts(pmids: list[str], use_cache: bool = True) -> dict[str, dict]:
    """Return {pmid: {title, abstract}} for the requested PMIDs.

    Cache-first: anything already in the local SQLite cache is reused.
    Missing PMIDs are fetched from NCBI E-utils in batches and persisted.
    """
    unique = list(dict.fromkeys(str(p).strip() for p in pmids if p))
    if not unique:
        return {}

    cached = _cache_get(unique) if use_cache else {}
    missing = [p for p in unique if p not in cached]
    if not missing:
        return cached

    fetched: dict[str, dict] = {}
    for i in range(0, len(missing), NCBI_BATCH):
        batch = missing[i : i + NCBI_BATCH]
        try:
            r = requests.get(
                NCBI_EFETCH,
                params={"db": "pubmed", "id": ",".join(batch),
                        "rettype": "abstract", "retmode": "xml"},
                timeout=30,
            )
            r.raise_for_status()
            fetched.update(_parse_pubmed_xml(r.text, batch))
        except Exception as e:
            print(f"[pubmed] fetch error for batch starting {batch[0]}: {e}", file=sys.stderr)
            for p in batch:
                fetched.setdefault(p, {"title": "", "abstract": ""})
        if i + NCBI_BATCH < len(missing):
            time.sleep(NCBI_SLEEP)

    _cache_put(fetched)
    cached.update(fetched)
    return cached


def cache_stats() -> dict:
    conn = _cache_conn()
    n = conn.execute("SELECT COUNT(*) FROM pubmed_abstract").fetchone()[0]
    conn.close()
    return {"cache_path": str(CACHE_PATH), "cached_abstracts": n}
