#!/usr/bin/env python3
"""
Deep-mine PubMed abstracts and GlyGen ptm_annotation for quantitative,
site-specific PTM changes between conditions. Extract:
  - Site + condition + direction (up/down) + fold change / % if mentioned
  - Sentences that co-mention a site AND a condition AND a direction
"""

import json, re, sys, time, requests

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ─── Regex patterns ───
SITE_PAT = re.compile(
    r'(?:'
    r'(?:Ser|Thr|Tyr|Lys|Asn|Cys|Arg|His|Gly|Ala|Val|Leu|Ile|Pro|Phe|Met|Trp|Asp|Glu|Gln|Trp)'
    r'[\s\-]?(\d{1,4})' 
    r'|'
    r'(?:p?[STYNKCRHDEGAQWFLMIPV])(\d{2,4})\b'
    r'|'
    r'(?:position|residue|site)\s+(\d{2,4})'
    r')',
    re.IGNORECASE
)

FOLD_PAT = re.compile(r'(\d+\.?\d*)\s*-?\s*fold', re.IGNORECASE)
PCT_PAT = re.compile(r'(\d+\.?\d*)\s*%', re.IGNORECASE)
RATIO_PAT = re.compile(r'ratio\s+(?:of\s+)?(\d+\.?\d*)', re.IGNORECASE)

UP_WORDS = [
    "increased", "upregulated", "up-regulated", "elevated", "enhanced",
    "higher", "more prevalent", "accumulation", "hyperphosphorylat",
    "overexpress", "gain", "promotes", "activates", "stimulates",
    "induces", "augment", "amplif", "enrich", "abundant",
]
DOWN_WORDS = [
    "decreased", "downregulated", "down-regulated", "reduced", "lower",
    "loss", "inhibits", "blocks", "prevents", "impaired", "deficient",
    "absent", "abolishes", "disrupts", "attenuate", "diminish",
    "deplet", "suppress",
]

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


def fetch_all_abstracts(pmids, batch_size=80):
    """Fetch PubMed abstracts and store full text per PMID."""
    out = {}
    pmids = list(set(pmids))
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]
        try:
            r = requests.get(NCBI_EFETCH, params={
                "db": "pubmed", "id": ",".join(batch),
                "rettype": "abstract", "retmode": "xml",
            }, timeout=30)
            if r.status_code != 200:
                continue
            xml = r.text
            for pmid in batch:
                blocks = re.findall(
                    r'<PMID[^>]*>' + re.escape(pmid) + r'</PMID>.*?</PubmedArticle>',
                    xml, re.DOTALL
                )
                if not blocks:
                    continue
                block = blocks[0]
                title_m = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', block, re.DOTALL)
                title = re.sub(r'<[^>]+>', '', title_m.group(1)) if title_m else ""
                abs_parts = re.findall(r'<AbstractText[^>]*>(.*?)</AbstractText>', block, re.DOTALL)
                abstract = " ".join(re.sub(r'<[^>]+>', '', p) for p in abs_parts)
                out[pmid] = {"title": title.strip(), "abstract": abstract.strip()}
        except Exception as e:
            print(f"  batch err: {e}", file=sys.stderr)
        time.sleep(0.35)
    return out


def extract_direction(text):
    tl = text.lower()
    for w in UP_WORDS:
        if w.lower() in tl:
            return "up"
    for w in DOWN_WORDS:
        if w.lower() in tl:
            return "down"
    return None


def extract_quant(text):
    """Extract fold-change or percentage from text."""
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


def match_terms(text, term_map):
    results = set()
    for pat, name in term_map.items():
        if re.search(pat, text, re.IGNORECASE):
            results.add(name)
    return sorted(results)


def mine_sentence(sent, protein_ac, protein_name, gene_name):
    """Mine a single sentence for site-specific quantitative PTM info."""
    # Must mention at least a site or a PTM
    sites = []
    for m in SITE_PAT.finditer(sent):
        pos = m.group(1) or m.group(2) or m.group(3)
        if pos:
            # Try to get the residue letter from context
            pre = sent[max(0, m.start()-4):m.start()]
            aa_m = re.search(r'(Ser|Thr|Tyr|Lys|Asn|Cys|Arg|Glu|Asp|Gln|His|[STYNKC])\s*$', pre, re.IGNORECASE)
            aa = aa_m.group(1) if aa_m else ""
            sites.append(f"{aa}{pos}".strip())

    diseases = match_terms(sent, DISEASE_MAP)
    tissues = match_terms(sent, TISSUE_MAP)
    ptms = match_terms(sent, PTM_MAP)
    direction = extract_direction(sent)
    quant = extract_quant(sent)

    # Only return if we have at least direction + (disease or tissue)
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


def main():
    # Load existing data
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/glygen_data.json") as f:
        gdata = json.load(f)

    # ─── PASS 1: Mine ptm_annotation free text (curated, highest quality) ───
    print("Pass 1: Mining curated ptm_annotation text...", file=sys.stderr)
    curated_hits = []
    for protein in gdata["proteins"]:
        ac = protein["uniprot_ac"]
        name = protein["display_name"]
        gene = protein.get("gene_name", "")
        for ann in protein.get("ptm_annotations_raw", []):
            sentences = re.split(r'(?<=[.!?])\s+', ann)
            for sent in sentences:
                if len(sent) < 20:
                    continue
                hit = mine_sentence(sent, ac, name, gene)
                if hit:
                    hit["source"] = "GlyGen curated (UniProtKB)"
                    hit["confidence"] = "high"
                    # Extract PMIDs from parent annotation
                    curated_hits.append(hit)

    print(f"  Found {len(curated_hits)} curated hits", file=sys.stderr)

    # ─── PASS 2: Re-fetch raw protein data and collect all evidence PMIDs ───
    print("Pass 2: Fetching per-site evidence PMIDs from GlyGen...", file=sys.stderr)

    PROTEINS = {p["uniprot_ac"]: p["display_name"] for p in gdata["proteins"]}
    all_pmids = set()
    site_to_pmids = {}  # (ac, site_pos) -> list of pmids

    for idx, (ac, name) in enumerate(PROTEINS.items()):
        print(f"  [{idx+1}/{len(PROTEINS)}] {ac}", file=sys.stderr)
        try:
            r = requests.get(f"https://api.glygen.org/protein/detail/{ac}/", timeout=60)
            if r.status_code != 200:
                continue
            raw = r.json()
        except:
            continue

        gene = ""
        for gn in raw.get("gene_names", []):
            if gn.get("type") == "recommended":
                gene = gn["name"]
                break

        for section in ["glycosylation", "phosphorylation"]:
            for entry in raw.get(section, []):
                pmids = [e["id"] for e in entry.get("evidence", []) if e.get("database") == "PubMed"]
                pos = entry.get("start_pos", entry.get("site_lbl", ""))
                key = (ac, str(pos), section)
                site_to_pmids[key] = pmids
                all_pmids.update(pmids)

        time.sleep(0.3)

    print(f"  Total unique PMIDs: {len(all_pmids)}", file=sys.stderr)

    # ─── PASS 3: Fetch and deep-mine PubMed abstracts ───
    print("Pass 3: Fetching PubMed abstracts...", file=sys.stderr)
    abstracts = fetch_all_abstracts(list(all_pmids))
    print(f"  Fetched {len(abstracts)} abstracts", file=sys.stderr)

    print("Pass 4: Mining abstracts for quantitative site-specific changes...", file=sys.stderr)
    abstract_hits = []

    for pmid, ab in abstracts.items():
        full_text = f"{ab['title']}. {ab['abstract']}"
        sentences = re.split(r'(?<=[.!?])\s+', full_text)

        # Find which protein(s) this PMID is evidence for
        linked_proteins = set()
        for (ac, pos, section), pmids in site_to_pmids.items():
            if pmid in pmids:
                linked_proteins.add(ac)

        for sent in sentences:
            if len(sent) < 25:
                continue
            for pac in linked_proteins:
                pname = PROTEINS.get(pac, pac)
                gene = pname.split("(")[0].strip() if "(" in pname else ""
                hit = mine_sentence(sent, pac, pname, gene)
                if hit:
                    hit["source"] = "PubMed abstract"
                    hit["confidence"] = "medium"
                    hit["pmid"] = pmid
                    hit["paper_title"] = ab["title"][:150]
                    abstract_hits.append(hit)

    print(f"  Found {len(abstract_hits)} abstract hits", file=sys.stderr)

    # ─── Deduplicate and merge ───
    all_hits = curated_hits + abstract_hits

    # Build summary
    conditions_count = {}
    for h in all_hits:
        for d in h["diseases"]:
            conditions_count[d] = conditions_count.get(d, 0) + 1

    sites_with_direction = [h for h in all_hits if h["direction"] in ("up", "down")]
    sites_with_quant = [h for h in all_hits if h.get("quantitative")]

    summary = {
        "total_hits": len(all_hits),
        "curated_hits": len(curated_hits),
        "abstract_hits": len(abstract_hits),
        "with_direction": len(sites_with_direction),
        "with_quantitative": len(sites_with_quant),
        "unique_conditions": len(conditions_count),
        "proteins_covered": len(set(h["protein_ac"] for h in all_hits)),
        "pmids_mined": len(abstracts),
    }

    output = {
        "summary": summary,
        "hits": all_hits,
    }

    outpath = "/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/quantitative.json"
    with open(outpath, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    import os
    print(f"\nDone. {os.path.getsize(outpath)} bytes", file=sys.stderr)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
