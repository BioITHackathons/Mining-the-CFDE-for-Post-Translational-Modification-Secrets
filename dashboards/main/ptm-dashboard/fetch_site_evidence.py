#!/usr/bin/env python3
"""
For each glycosylation and phosphorylation site in our protein set:
1. Collect all PubMed IDs linked as evidence
2. Fetch abstracts from NCBI E-utilities
3. NLP-extract tissue, condition/disease, and experimental context from each abstract
4. Produce a unified site-level table with condition + tissue annotations
"""

import json, re, sys, time, requests
from collections import defaultdict

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Disease vocabulary for abstract mining
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


def fetch_abstracts_batch(pmids, batch_size=50):
    """Fetch PubMed abstracts in batches via NCBI E-utilities."""
    abstracts = {}
    pmid_list = list(set(pmids))
    for i in range(0, len(pmid_list), batch_size):
        batch = pmid_list[i:i+batch_size]
        ids = ",".join(batch)
        try:
            r = requests.get(NCBI_EFETCH, params={
                "db": "pubmed",
                "id": ids,
                "rettype": "abstract",
                "retmode": "xml",
            }, timeout=30)
            if r.status_code == 200:
                # Parse XML for abstracts and titles
                text = r.text
                # Simple XML extraction
                for pmid in batch:
                    # Find the article block for this PMID
                    pat = re.compile(
                        r'<PMID[^>]*>' + re.escape(pmid) + r'</PMID>.*?'
                        r'(?:<ArticleTitle>(.*?)</ArticleTitle>)?.*?'
                        r'(?:<AbstractText[^>]*>(.*?)</AbstractText>)?',
                        re.DOTALL
                    )
                    m = pat.search(text)
                    title = ""
                    abstract = ""
                    if m:
                        title = re.sub(r'<[^>]+>', '', m.group(1) or "")
                        abstract = re.sub(r'<[^>]+>', '', m.group(2) or "")
                    # Sometimes there are multiple AbstractText blocks
                    all_abs = re.findall(
                        r'<PMID[^>]*>' + re.escape(pmid) + r'</PMID>.*?</PubmedArticle>',
                        text, re.DOTALL
                    )
                    if all_abs:
                        abs_texts = re.findall(r'<AbstractText[^>]*>(.*?)</AbstractText>', all_abs[0], re.DOTALL)
                        if abs_texts:
                            abstract = " ".join(re.sub(r'<[^>]+>', '', t) for t in abs_texts)
                        title_m = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', all_abs[0], re.DOTALL)
                        if title_m:
                            title = re.sub(r'<[^>]+>', '', title_m.group(1))
                    abstracts[pmid] = {"title": title.strip(), "abstract": abstract.strip()}
        except Exception as e:
            print(f"  Error fetching batch starting {batch[0]}: {e}", file=sys.stderr)
        time.sleep(0.4)  # NCBI rate limit
    return abstracts


def mine_text(text):
    """Extract diseases, tissues, and experimental methods from text."""
    diseases = set()
    tissues = set()
    methods = set()
    text_lower = text.lower() if text else ""

    for pat, name in DISEASE_TERMS.items():
        if re.search(pat, text, re.IGNORECASE):
            diseases.add(name)
    for pat, name in TISSUE_TERMS.items():
        if re.search(pat, text, re.IGNORECASE):
            tissues.add(name)
    for pat, name in EXPERIMENTAL_TERMS.items():
        if re.search(pat, text, re.IGNORECASE):
            methods.add(name)

    return sorted(diseases), sorted(tissues), sorted(methods)


def main():
    # Load existing GlyGen data
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/glygen_data.json") as f:
        data = json.load(f)

    # Collect all unique PMIDs across all sites
    all_pmids = set()
    site_pmid_map = []  # (protein_ac, protein_name, site_lbl, ptm_type, glycan, pmid)

    for protein in data["proteins"]:
        ac = protein["uniprot_ac"]
        name = protein["display_name"]

        # Glycosylation sites
        for site in protein.get("glycosylation_sites", []):
            # We need full evidence - re-fetch from the raw data
            pass

        # Phosphorylation sites
        for site in protein.get("phosphorylation_sites", []):
            pass

    # Actually we need the full raw data with evidence PMIDs per site.
    # Let's re-read the full glygen_data.json which has evidence_count but not the actual PMIDs.
    # We need to go back to the raw API data. Let's fetch site details for key proteins.

    # Focused approach: For the proteins with the most annotations, fetch site-level detail
    # from GlyGen's site detail API to get per-site PMIDs

    # First, collect all PMIDs from publications, ptm_annotations, etc.
    print("Collecting PMIDs from protein-level data...", file=sys.stderr)

    # We'll use the full raw API data - fetch protein detail again for the top proteins
    # and extract per-site evidence PMIDs
    TOP_PROTEINS = [
        "P05067", "P10636", "P37840", "P02649", "Q9NZC2",
        "P56817", "P49768", "P04156", "P21802", "P10909",
        "P07602", "P02751", "P08473", "Q92673", "P20933",
        "P78509", "Q99700", "P49810", "P51693", "Q06481",
        "P06396", "P22301", "P09601", "P08254", "P16422",
        "P21589", "Q13421", "P01130", "P15169", "P22891",
    ]

    site_records = []

    for idx, ac in enumerate(TOP_PROTEINS):
        pname = {p["uniprot_ac"]: p["display_name"] for p in data["proteins"]}.get(ac, ac)
        print(f"[{idx+1}/{len(TOP_PROTEINS)}] Fetching raw detail for {ac} ({pname})...", file=sys.stderr)

        try:
            r = requests.get(f"https://api.glygen.org/protein/detail/{ac}/", timeout=60)
            if r.status_code != 200:
                continue
            raw = r.json()
        except Exception as e:
            print(f"  Error: {e}", file=sys.stderr)
            continue

        gene_name = ""
        for gn in raw.get("gene_names", []):
            if gn.get("type") == "recommended":
                gene_name = gn["name"]
                break

        protein_name = ""
        for pn in raw.get("protein_names", []):
            if pn.get("type") == "recommended":
                protein_name = pn["name"]
                break

        # Glycosylation sites with per-site evidence
        for g in raw.get("glycosylation", []):
            pmids = [e["id"] for e in g.get("evidence", []) if e.get("database") == "PubMed"]
            all_pmids.update(pmids)
            site_records.append({
                "uniprot_ac": ac,
                "display_name": pname,
                "gene_name": gene_name,
                "protein_name": protein_name,
                "site_lbl": g.get("site_lbl", f"pos{g.get('start_pos', '?')}"),
                "start_pos": g.get("start_pos", ""),
                "residue": g.get("residue", g.get("site_seq", "")),
                "ptm_type": "Glycosylation",
                "ptm_subtype": g.get("type", ""),  # N-linked, O-linked
                "glycan_subtype": g.get("subtype", ""),
                "glytoucan_ac": g.get("glytoucan_ac", ""),
                "category": g.get("site_category", ""),
                "pmids": pmids,
            })

        # Phosphorylation sites with evidence
        for p in raw.get("phosphorylation", []):
            pmids = [e["id"] for e in p.get("evidence", []) if e.get("database") == "PubMed"]
            all_pmids.update(pmids)
            site_records.append({
                "uniprot_ac": ac,
                "display_name": pname,
                "gene_name": gene_name,
                "protein_name": protein_name,
                "site_lbl": p.get("site_lbl", f"pos{p.get('start_pos', '?')}"),
                "start_pos": p.get("start_pos", ""),
                "residue": p.get("residue", ""),
                "ptm_type": "Phosphorylation",
                "ptm_subtype": p.get("comment", ""),
                "glycan_subtype": "",
                "glytoucan_ac": "",
                "category": "reported",
                "pmids": pmids,
            })

        time.sleep(0.5)

    print(f"\nTotal site records: {len(site_records)}", file=sys.stderr)
    print(f"Total unique PMIDs to fetch: {len(all_pmids)}", file=sys.stderr)

    # Fetch abstracts
    print("Fetching PubMed abstracts...", file=sys.stderr)
    abstracts = fetch_abstracts_batch(list(all_pmids))
    print(f"Fetched {len(abstracts)} abstracts", file=sys.stderr)

    # Mine each abstract for disease/tissue/method
    print("Mining abstracts for disease/tissue/method context...", file=sys.stderr)
    abstract_annotations = {}
    for pmid, ab in abstracts.items():
        combined = f"{ab['title']} {ab['abstract']}"
        diseases, tissues, methods = mine_text(combined)
        abstract_annotations[pmid] = {
            "title": ab["title"],
            "diseases": diseases,
            "tissues": tissues,
            "methods": methods,
        }

    # Build final site-level records with condition/tissue from evidence
    print("Building site-level condition/tissue annotations...", file=sys.stderr)
    enriched_sites = []

    for rec in site_records:
        # Aggregate diseases/tissues from all linked abstracts
        all_diseases = set()
        all_tissues = set()
        all_methods = set()
        evidence_details = []

        for pmid in rec["pmids"]:
            ann = abstract_annotations.get(pmid, {})
            if ann:
                all_diseases.update(ann.get("diseases", []))
                all_tissues.update(ann.get("tissues", []))
                all_methods.update(ann.get("methods", []))
                evidence_details.append({
                    "pmid": pmid,
                    "title": ann.get("title", ""),
                    "diseases": ann.get("diseases", []),
                    "tissues": ann.get("tissues", []),
                    "methods": ann.get("methods", []),
                })

        enriched_sites.append({
            **rec,
            "conditions": sorted(all_diseases),
            "tissues": sorted(all_tissues),
            "methods": sorted(all_methods),
            "evidence": evidence_details,
        })

    # Summary stats
    sites_with_conditions = sum(1 for s in enriched_sites if s["conditions"])
    sites_with_tissues = sum(1 for s in enriched_sites if s["tissues"])
    unique_conditions = set()
    unique_tissues = set()
    for s in enriched_sites:
        unique_conditions.update(s["conditions"])
        unique_tissues.update(s["tissues"])

    summary = {
        "total_site_records": len(enriched_sites),
        "sites_with_condition_context": sites_with_conditions,
        "sites_with_tissue_context": sites_with_tissues,
        "unique_conditions": len(unique_conditions),
        "unique_tissues": len(unique_tissues),
        "unique_pmids_mined": len(abstract_annotations),
        "proteins_covered": len(set(s["uniprot_ac"] for s in enriched_sites)),
    }

    output = {
        "summary": summary,
        "sites": enriched_sites,
    }

    outpath = "/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/site_evidence.json"
    with open(outpath, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    print(f"\nDone. Wrote {outpath}", file=sys.stderr)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
