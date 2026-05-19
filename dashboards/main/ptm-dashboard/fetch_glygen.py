#!/usr/bin/env python3
"""
Fetch GlyGen protein detail data for neurodegenerative-disease-relevant glycoproteins,
extract structured + free-text PTM/disease/tissue annotations, and produce a unified JSON
for the interactive dashboard.
"""

import json, re, sys, time, requests

# Target proteins: neurodegeneration-relevant glycoproteins + key PTM substrates
PROTEINS = {
    "P05067": "APP (Amyloid Precursor Protein)",
    "P10636": "MAPT (Tau)",
    "P37840": "SNCA (Alpha-synuclein)",
    "P02649": "APOE (Apolipoprotein E)",
    "Q9NZC2": "TREM2",
    "P01130": "LDLR (LDL Receptor)",
    "P56817": "BACE1 (Beta-secretase 1)",
    "P49768": "PSEN1 (Presenilin-1)",
    "P49810": "PSEN2 (Presenilin-2)",
    "P04156": "PRNP (Prion Protein)",
    "P09601": "HMOX1 (Heme Oxygenase 1)",
    "P78509": "RELN (Reelin)",
    "P21802": "FGFR2",
    "Q99700": "ATXN2 (Ataxin-2)",
    "P51693": "APLP1 (APP-like Protein 1)",
    "Q06481": "APLP2 (APP-like Protein 2)",
    "P20933": "NPC1 (Niemann-Pick C1)",
    "Q92673": "SORL1 (Sortilin-related receptor)",
    "P08254": "MMP3 (Matrix Metalloproteinase 3)",
    "P16422": "EPCAM",
    "P22301": "IL10 (Interleukin-10)",
    "P10909": "CLU (Clusterin)",
    "P07602": "PSAP (Prosaposin)",
    "P02751": "FN1 (Fibronectin)",
    "P06396": "GSN (Gelsolin)",
    "Q13421": "MSLN (Mesothelin)",
    "P21589": "NT5E (CD73)",
    "P08473": "MME (Neprilysin)",
    "P15169": "CPN1 (Carboxypeptidase N)",
    "P22891": "PROZ (Protein Z)",
}

def fetch_protein(uniprot_ac):
    """Fetch full protein detail from GlyGen API."""
    url = f"https://api.glygen.org/protein/detail/{uniprot_ac}/"
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                return r.json()
            print(f"  Warning: {uniprot_ac} returned status {r.status_code}", file=sys.stderr)
        except Exception as e:
            print(f"  Retry {attempt+1} for {uniprot_ac}: {e}", file=sys.stderr)
            time.sleep(2)
    return None

def extract_ptm_annotations(annotations):
    """Parse ptm_annotation free-text into structured disease/site/tissue mentions."""
    results = []
    if not annotations:
        return results

    # Patterns for extracting site mentions
    site_pat = re.compile(r'(?:Ser|Thr|Tyr|Lys|Asn|Cys|Arg|His|Gly|Ala|Val|Leu|Ile|Pro|Phe|Met|Trp|Asp|Glu|Gln)[\-\s]?(\d+)', re.IGNORECASE)
    site_pat2 = re.compile(r'([STYNKCRHDEGAQWFLMIPV])(\d{2,4})\b')

    # Disease keywords
    disease_keywords = {
        "alzheimer": "Alzheimer's disease",
        "parkinson": "Parkinson's disease",
        "huntington": "Huntington's disease",
        "amyotrophic lateral sclerosis": "ALS",
        "ALS": "ALS",
        "prion": "Prion disease",
        "dementia": "Dementia",
        "neurodegenerat": "Neurodegenerative disease",
        "cancer": "Cancer",
        "tumor": "Cancer",
        "carcinoma": "Cancer",
        "diabetes": "Diabetes",
        "inflammatory": "Inflammatory disease",
        "autoimmune": "Autoimmune disease",
        "cardiovascular": "Cardiovascular disease",
        "atherosclerosis": "Atherosclerosis",
        "fibrosis": "Fibrosis",
        "multiple sclerosis": "Multiple sclerosis",
        "epilepsy": "Epilepsy",
        "schizophrenia": "Schizophrenia",
        "stroke": "Stroke",
        "ischemia": "Ischemic disease",
        "amyloid": "Amyloidosis",
        "tauopathy": "Tauopathy",
        "lewy bod": "Lewy body disease",
        "frontotemporal": "Frontotemporal dementia",
        "cerebral amyloid angiopathy": "Cerebral amyloid angiopathy",
        "down syndrome": "Down syndrome",
    }

    # Tissue/cell type keywords
    tissue_keywords = {
        "neuron": "Neuron",
        "neuronal": "Neuron",
        "brain": "Brain",
        "cerebr": "Brain (cerebral)",
        "hippocamp": "Hippocampus",
        "cortex": "Cortex",
        "cortical": "Cortex",
        "astrocyte": "Astrocyte",
        "microglia": "Microglia",
        "oligodendrocyte": "Oligodendrocyte",
        "synap": "Synapse",
        "axon": "Axon",
        "dendrit": "Dendrite",
        "dopaminer": "Dopaminergic neuron",
        "liver": "Liver",
        "hepat": "Liver",
        "kidney": "Kidney",
        "renal": "Kidney",
        "plasma": "Blood plasma",
        "serum": "Blood serum",
        "erythrocyte": "Erythrocyte",
        "platelet": "Platelet",
        "endotheli": "Endothelium",
        "epitheli": "Epithelium",
        "mitochondri": "Mitochondria",
        "lysosom": "Lysosome",
        "endoplasmic reticulum": "Endoplasmic reticulum",
        "golgi": "Golgi apparatus",
        "membrane": "Cell membrane",
        "extracellular": "Extracellular",
        "cerebrospinal fluid": "CSF",
        "CSF": "CSF",
        "substantia nigra": "Substantia nigra",
        "striatum": "Striatum",
        "spinal cord": "Spinal cord",
    }

    # PTM type keywords
    ptm_keywords = {
        "phosphorylat": "Phosphorylation",
        "glycosylat": "Glycosylation",
        "n-glyco": "N-glycosylation",
        "o-glyco": "O-glycosylation",
        "ubiquitin": "Ubiquitination",
        "sumoylat": "SUMOylation",
        "acetylat": "Acetylation",
        "methylat": "Methylation",
        "palmitoyl": "Palmitoylation",
        "sulfat": "Sulfation",
        "nitrosylat": "S-nitrosylation",
        "sialylat": "Sialylation",
        "fucosylat": "Fucosylation",
        "mannosylat": "Mannosylation",
        "proteolytic": "Proteolytic cleavage",
        "cleavage": "Proteolytic cleavage",
        "truncat": "Truncation",
        "aggregat": "Aggregation",
        "glycation": "Glycation",
    }

    # Direction keywords
    up_words = ["increased", "upregulated", "up-regulated", "elevated", "enhanced", "higher",
                "more prevalent", "accumulation", "hyperphosphorylat", "overexpress", "gain",
                "promotes", "activates", "stimulates", "induces"]
    down_words = ["decreased", "downregulated", "down-regulated", "reduced", "lower", "loss",
                  "inhibits", "blocks", "prevents", "impaired", "deficient", "absent",
                  "abolishes", "disrupts"]

    for ann in annotations:
        text = ann.get("annotation", "")
        pmids = [e["id"] for e in ann.get("evidence", []) if e.get("database") == "PubMed"]

        # Split into sentences for finer-grained extraction
        sentences = re.split(r'(?<=[.!?])\s+', text)

        for sent in sentences:
            sent_lower = sent.lower()

            # Extract sites
            sites = []
            for m in site_pat.finditer(sent):
                aa = m.group(0).split("-")[0].split()[0] if "-" in m.group(0) or " " in m.group(0) else m.group(0)[:3]
                pos = m.group(1)
                sites.append(f"{aa}-{pos}")
            for m in site_pat2.finditer(sent):
                sites.append(f"{m.group(1)}{m.group(2)}")

            # Extract diseases
            diseases_found = []
            for kw, disease_name in disease_keywords.items():
                if kw.lower() in sent_lower:
                    diseases_found.append(disease_name)

            # Extract tissues
            tissues_found = []
            for kw, tissue_name in tissue_keywords.items():
                if kw.lower() in sent_lower:
                    tissues_found.append(tissue_name)

            # Extract PTM types
            ptms_found = []
            for kw, ptm_name in ptm_keywords.items():
                if kw.lower() in sent_lower:
                    ptms_found.append(ptm_name)

            # Extract direction
            direction = "unknown"
            for w in up_words:
                if w.lower() in sent_lower:
                    direction = "up"
                    break
            if direction == "unknown":
                for w in down_words:
                    if w.lower() in sent_lower:
                        direction = "down"
                        break

            # Only keep sentences that have at least one PTM or site mention
            if ptms_found or sites:
                results.append({
                    "sentence": sent.strip(),
                    "sites": list(set(sites)),
                    "diseases": list(set(diseases_found)),
                    "tissues": list(set(tissues_found)),
                    "ptm_types": list(set(ptms_found)),
                    "direction": direction,
                    "pmids": pmids,
                })

    return results

def process_protein(uniprot_ac, display_name, raw):
    """Process a single protein's GlyGen data into the dashboard format."""
    protein_data = {
        "uniprot_ac": uniprot_ac,
        "display_name": display_name,
        "gene_name": raw.get("uniprot", {}).get("uniprot_id", ""),
        "protein_name": "",
        "length": raw.get("sequence", {}).get("length", 0),
        "species": "",
    }

    # Protein name
    names = raw.get("protein_names", [])
    for n in names:
        if n.get("type") == "recommended":
            protein_data["protein_name"] = n["name"]
            break

    # Species
    species = raw.get("species", [])
    if species:
        protein_data["species"] = species[0].get("name", "")

    # Structured glycosylation sites
    glyco_sites = []
    for g in raw.get("glycosylation", []):
        glyco_sites.append({
            "start_pos": g.get("start_pos", g.get("start_aa", "")),
            "end_pos": g.get("end_pos", g.get("end_aa", "")),
            "residue": g.get("residue", g.get("site_seq", "")),
            "site_lbl": g.get("site_lbl", ""),
            "type": g.get("type", ""),
            "subtype": g.get("subtype", ""),
            "glytoucan_ac": g.get("glytoucan_ac", ""),
            "category": g.get("site_category", ""),
            "evidence_count": len(g.get("evidence", [])),
        })
    protein_data["glycosylation_sites"] = glyco_sites
    protein_data["glycosylation_count"] = len(glyco_sites)

    # Structured phosphorylation sites
    phospho_sites = []
    for p in raw.get("phosphorylation", []):
        phospho_sites.append({
            "start_pos": p.get("start_pos", ""),
            "residue": p.get("residue", ""),
            "site_lbl": p.get("site_lbl", ""),
            "comment": p.get("comment", ""),
            "evidence_count": len(p.get("evidence", [])),
        })
    protein_data["phosphorylation_sites"] = phospho_sites
    protein_data["phosphorylation_count"] = len(phospho_sites)

    # Disease associations (protein-level)
    diseases = []
    for d in raw.get("disease", []):
        rn = d.get("recommended_name", {})
        diseases.append({
            "disease_id": d.get("disease_id", ""),
            "name": rn.get("name", ""),
            "description": rn.get("description", ""),
        })
    protein_data["diseases"] = diseases

    # Expression in disease
    expr_disease = []
    for ed in raw.get("expression_disease", []):
        expr_disease.append({
            "disease_name": ed.get("name", ""),
            "trend": ed.get("trend", ""),
            "significant": ed.get("significant", ""),
        })
    protein_data["expression_disease"] = expr_disease

    # Expression in tissue
    expr_tissue = []
    for et in raw.get("expression_tissue", []):
        t = et.get("tissue", {})
        expr_tissue.append({
            "tissue_name": t.get("name", ""),
            "tissue_id": f"{t.get('namespace', '')}:{t.get('id', '')}",
            "score": float(et.get("score", 0)),
            "present": et.get("present", ""),
        })
    protein_data["expression_tissue"] = expr_tissue

    # SNVs with disease links and glycoeffects
    snv_disease = []
    for s in raw.get("snv", []):
        if s.get("disease") or s.get("glycoeffect"):
            snv_entry = {
                "site_lbl": s.get("site_lbl", ""),
                "start_pos": s.get("start_pos", ""),
                "sequence_org": s.get("sequence_org", ""),
                "sequence_mut": s.get("sequence_mut", ""),
                "glycoeffect": s.get("glycoeffect", []),
                "comment": s.get("comment", ""),
                "keywords": s.get("keywords", []),
                "diseases": [],
            }
            for sd in s.get("disease", []):
                rn = sd.get("recommended_name", {})
                snv_entry["diseases"].append({
                    "disease_id": sd.get("disease_id", ""),
                    "name": rn.get("name", ""),
                })
            snv_disease.append(snv_entry)
    protein_data["snv_with_disease_or_glycoeffect"] = snv_disease

    # PTM annotations (free text) - NLP extraction
    ptm_annotations_raw = raw.get("ptm_annotation", [])
    protein_data["ptm_annotations_raw"] = [a.get("annotation", "") for a in ptm_annotations_raw]
    protein_data["ptm_annotations_extracted"] = extract_ptm_annotations(ptm_annotations_raw)

    # Publications count
    protein_data["publication_count"] = len(raw.get("publication", []))

    return protein_data

def main():
    all_data = []
    total = len(PROTEINS)
    for i, (ac, name) in enumerate(PROTEINS.items()):
        print(f"[{i+1}/{total}] Fetching {ac} ({name})...", file=sys.stderr)
        raw = fetch_protein(ac)
        if raw:
            processed = process_protein(ac, name, raw)
            all_data.append(processed)
            print(f"  -> {processed['glycosylation_count']} glyco sites, {processed['phosphorylation_count']} phospho sites, {len(processed['diseases'])} diseases, {len(processed['ptm_annotations_extracted'])} extracted annotations", file=sys.stderr)
        else:
            print(f"  -> FAILED to fetch {ac}", file=sys.stderr)
        time.sleep(0.5)  # Be polite to the API

    # Summary stats
    summary = {
        "total_proteins": len(all_data),
        "total_glyco_sites": sum(p["glycosylation_count"] for p in all_data),
        "total_phospho_sites": sum(p["phosphorylation_count"] for p in all_data),
        "total_diseases": len(set(d["disease_id"] for p in all_data for d in p["diseases"])),
        "total_extracted_annotations": sum(len(p["ptm_annotations_extracted"]) for p in all_data),
        "total_snv_disease_links": sum(len(p["snv_with_disease_or_glycoeffect"]) for p in all_data),
        "data_source": "GlyGen API (https://api.glygen.org/)",
        "cfde_integration": "GlyGen is a CFDE DCC providing glycosylation data via C2M2",
    }

    output = {
        "summary": summary,
        "proteins": all_data,
    }

    outpath = "/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/glygen_data.json"
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nDone. Wrote {outpath}", file=sys.stderr)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
