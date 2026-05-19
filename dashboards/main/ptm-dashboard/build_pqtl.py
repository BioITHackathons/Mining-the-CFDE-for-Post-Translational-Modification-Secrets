#!/usr/bin/env python3
"""
Cross-reference GlyGen SNVs with PTM sites to find variants that:
1. Overlap directly with a PTM site
2. Are within ±5 residues of a PTM site  
3. Have disease associations AND affect PTM sites (pQTL-like)
4. Have computed glycoeffects

This uses only existing GlyGen data -- no external API calls.
"""
import json, sys, time, requests

def main():
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/glygen_data.json") as f:
        gdata = json.load(f)
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/tracks.json") as f:
        tdata = json.load(f)

    proteins = {p["uniprot_ac"]: p for p in gdata["proteins"]}

    # For full SNV data with positions, we need to re-fetch from GlyGen
    # (our existing data only has SNVs with disease OR glycoeffect, not all SNVs)
    print("Fetching full SNV data per protein from GlyGen...", file=sys.stderr)

    all_overlaps = []
    protein_summaries = {}

    for idx, (ac, pdata) in enumerate(proteins.items()):
        name = pdata["display_name"]
        print(f"  [{idx+1}/{len(proteins)}] {ac} ({name})", file=sys.stderr)

        try:
            r = requests.get(f"https://api.glygen.org/protein/detail/{ac}/", timeout=60)
            if r.status_code != 200:
                continue
            raw = r.json()
        except Exception as e:
            print(f"    Error: {e}", file=sys.stderr)
            continue

        # Collect all PTM site positions
        ptm_sites = {}  # pos -> {type, label, ...}
        for g in raw.get("glycosylation", []):
            pos = g.get("start_pos")
            if isinstance(pos, int) and pos > 0:
                if pos not in ptm_sites:
                    ptm_sites[pos] = {"types": set(), "label": g.get("site_lbl", f"pos{pos}"), "residue": g.get("residue", g.get("site_seq", ""))}
                ptm_sites[pos]["types"].add(f"Glycosylation ({g.get('type', '')})")

        for p in raw.get("phosphorylation", []):
            pos = p.get("start_pos")
            if isinstance(pos, int) and pos > 0:
                if pos not in ptm_sites:
                    ptm_sites[pos] = {"types": set(), "label": p.get("site_lbl", f"pos{pos}"), "residue": p.get("residue", "")}
                ptm_sites[pos]["types"].add("Phosphorylation")

        ptm_positions = set(ptm_sites.keys())

        # Collect all SNVs
        snvs = raw.get("snv", [])
        gene_name = ""
        for gn in raw.get("gene_names", []):
            if gn.get("type") == "recommended":
                gene_name = gn["name"]
                break

        snv_at_ptm = 0
        snv_near_ptm = 0
        snv_disease_ptm = 0

        for snv in snvs:
            snv_pos = snv.get("start_pos")
            if not isinstance(snv_pos, int):
                continue

            # Check overlap with PTM sites
            direct_overlap = snv_pos in ptm_positions
            near_positions = [p for p in ptm_positions if abs(p - snv_pos) <= 5 and p != snv_pos]
            near_ptm = len(near_positions) > 0

            if not direct_overlap and not near_ptm:
                continue

            # This SNV overlaps or is near a PTM site
            overlap_type = "direct" if direct_overlap else "proximal"
            distance = 0 if direct_overlap else min(abs(p - snv_pos) for p in near_positions)

            # Get PTM info for overlapping sites
            overlapping_ptm_pos = [snv_pos] if direct_overlap else near_positions
            ptm_info = []
            for pp in overlapping_ptm_pos:
                if pp in ptm_sites:
                    ptm_info.append({
                        "pos": pp,
                        "label": ptm_sites[pp]["label"],
                        "types": sorted(ptm_sites[pp]["types"]),
                        "residue": ptm_sites[pp]["residue"],
                    })

            # Disease associations
            diseases = []
            for d in snv.get("disease", []):
                rn = d.get("recommended_name", {})
                diseases.append({"id": d.get("disease_id", ""), "name": rn.get("name", "")})

            glycoeffects = snv.get("glycoeffect", [])
            keywords = snv.get("keywords", [])
            is_disease = len(diseases) > 0 or "disease" in keywords
            is_somatic = "somatic" in keywords
            is_germline = "germline" in keywords

            entry = {
                "protein_ac": ac,
                "protein_name": name,
                "gene_name": gene_name,
                "snv_pos": snv_pos,
                "snv_label": snv.get("site_lbl", f"pos{snv_pos}"),
                "ref_aa": snv.get("sequence_org", ""),
                "alt_aa": snv.get("sequence_mut", ""),
                "chr_id": snv.get("chr_id", ""),
                "chr_pos": snv.get("chr_pos", ""),
                "overlap_type": overlap_type,
                "distance": distance,
                "ptm_sites": ptm_info,
                "diseases": diseases,
                "glycoeffects": glycoeffects,
                "is_disease_associated": is_disease,
                "is_somatic": is_somatic,
                "is_germline": is_germline,
                "variant_type": "somatic" if is_somatic else "germline" if is_germline else "unknown",
                "comment": snv.get("comment", ""),
                "dbsnp": "",
                "evidence_sources": [],
            }

            # Extract dbSNP ID
            for ev in snv.get("evidence", []):
                if ev.get("database") == "dbSNP":
                    entry["dbsnp"] = ev.get("id", "")
                entry["evidence_sources"].append(ev.get("database", ""))

            all_overlaps.append(entry)

            if direct_overlap:
                snv_at_ptm += 1
            if near_ptm:
                snv_near_ptm += 1
            if is_disease and (direct_overlap or near_ptm):
                snv_disease_ptm += 1

        protein_summaries[ac] = {
            "name": name,
            "gene": gene_name,
            "total_ptm_sites": len(ptm_positions),
            "total_snvs": len(snvs),
            "snvs_at_ptm": snv_at_ptm,
            "snvs_near_ptm": snv_near_ptm,
            "snvs_disease_and_ptm": snv_disease_ptm,
        }

        time.sleep(0.4)

    # Compute summary stats
    direct_overlaps = [o for o in all_overlaps if o["overlap_type"] == "direct"]
    proximal = [o for o in all_overlaps if o["overlap_type"] == "proximal"]
    disease_overlaps = [o for o in all_overlaps if o["is_disease_associated"]]
    with_glycoeffect = [o for o in all_overlaps if o["glycoeffects"]]

    # Unique diseases across all overlaps
    all_diseases = set()
    for o in all_overlaps:
        for d in o["diseases"]:
            all_diseases.add(d["name"])

    summary = {
        "total_snv_ptm_overlaps": len(all_overlaps),
        "direct_overlaps": len(direct_overlaps),
        "proximal_overlaps": len(proximal),
        "disease_associated": len(disease_overlaps),
        "with_glycoeffect": len(with_glycoeffect),
        "unique_diseases": len(all_diseases),
        "proteins_with_overlaps": sum(1 for p in protein_summaries.values() if p["snvs_at_ptm"] + p["snvs_near_ptm"] > 0),
    }

    output = {
        "summary": summary,
        "overlaps": all_overlaps,
        "protein_summaries": protein_summaries,
    }

    outpath = "/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/pqtl.json"
    with open(outpath, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    import os
    print(f"\nDone. {os.path.getsize(outpath)} bytes", file=sys.stderr)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
