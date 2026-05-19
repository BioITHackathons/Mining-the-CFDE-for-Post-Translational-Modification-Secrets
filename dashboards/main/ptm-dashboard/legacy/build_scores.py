#!/usr/bin/env python3
"""
Build evidence-weighted logFC scores per protein-site per condition.

Scoring logic:
- Each "up" hit from a PubMed abstract = +1
- Each "down" hit from a PubMed abstract = -1
- Curated hit (from UniProtKB annotation) = +/- 2 (higher confidence weight)
- If explicit fold-change is in the abstract, use that as the score for that hit
- Aggregate: sum scores across all hits → net evidence score (pseudo-logFC)
- Also compute: total evidence count, # papers up, # papers down, max fold-change if any

This gives every site a quantitative score per condition that can be plotted
as logFC-style bars in the track viewer.
"""
import json, re, sys

def main():
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/quantitative.json") as f:
        qdata = json.load(f)
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/tracks.json") as f:
        tdata = json.load(f)

    # Build score index: (protein_ac, site_pos, condition) -> score data
    # We need to match quantitative hits to specific site positions

    # First pass: aggregate directional hits by protein + condition
    # (many hits don't have specific site positions, so we also track protein-level)
    protein_cond_scores = {}  # (ac, condition) -> {up: n, down: n, ...}
    site_cond_scores = {}     # (ac, pos, condition) -> {score, up, down, evidence, ...}

    for hit in qdata["hits"]:
        ac = hit["protein_ac"]
        direction = hit["direction"]
        confidence = hit.get("confidence", "medium")
        quant = hit.get("quantitative")
        pmid = hit.get("pmid", "")
        sentence = hit.get("sentence", "")[:150]

        # Score per hit
        if quant and quant.get("value"):
            # Use the actual fold-change value
            val = quant["value"]
            if quant["type"] == "percent":
                val = val / 100  # Normalize percentages
            hit_score = val if direction == "up" else -val
        elif confidence == "high":
            hit_score = 2.0 if direction == "up" else -2.0
        else:
            hit_score = 1.0 if direction == "up" else -1.0

        for condition in hit.get("diseases", []):
            # Skip generic conditions for cleaner comparisons
            if condition in ("Control", "Normal/healthy"):
                continue

            # Protein-level aggregation
            pk = (ac, condition)
            if pk not in protein_cond_scores:
                protein_cond_scores[pk] = {"up": 0, "down": 0, "score": 0, "evidence": []}
            pcs = protein_cond_scores[pk]
            if direction == "up":
                pcs["up"] += 1
            else:
                pcs["down"] += 1
            pcs["score"] += hit_score
            if len(pcs["evidence"]) < 3:
                pcs["evidence"].append({"pmid": pmid, "sentence": sentence, "score": hit_score})

            # Site-level aggregation (if specific sites are mentioned)
            for site_str in hit.get("sites", []):
                m = re.search(r'(\d+)', site_str)
                if not m:
                    continue
                pos = int(m.group(1))
                sk = (ac, pos, condition)
                if sk not in site_cond_scores:
                    site_cond_scores[sk] = {
                        "up": 0, "down": 0, "score": 0.0,
                        "evidence": [], "max_fc": None, "site_str": site_str,
                    }
                scs = site_cond_scores[sk]
                if direction == "up":
                    scs["up"] += 1
                else:
                    scs["down"] += 1
                scs["score"] += hit_score
                if quant and quant.get("value"):
                    if scs["max_fc"] is None or abs(quant["value"]) > abs(scs["max_fc"]):
                        scs["max_fc"] = quant["value"] if direction == "up" else -quant["value"]
                if len(scs["evidence"]) < 3:
                    scs["evidence"].append({"pmid": pmid, "sentence": sentence, "score": hit_score})

    # Now merge scores into the tracks data
    # For each protein → site, add condition scores
    for ac, pdata in tdata["proteins"].items():
        for site in pdata["sites"]:
            pos = site["pos"]
            site["scores"] = {}

            for condition in site.get("conditions", {}).keys():
                sk = (ac, pos, condition)
                if sk in site_cond_scores:
                    scs = site_cond_scores[sk]
                    site["scores"][condition] = {
                        "score": round(scs["score"], 2),
                        "up": scs["up"],
                        "down": scs["down"],
                        "max_fc": round(scs["max_fc"], 2) if scs["max_fc"] is not None else None,
                        "evidence": scs["evidence"],
                    }
                else:
                    # No site-specific quantitative hit; use protein-level as fallback
                    pk = (ac, condition)
                    if pk in protein_cond_scores:
                        pcs = protein_cond_scores[pk]
                        # Use a diluted protein-level score (divide by number of sites)
                        n_sites = max(len(pdata["sites"]), 1)
                        site["scores"][condition] = {
                            "score": round(pcs["score"] / n_sites, 2),
                            "up": pcs["up"],
                            "down": pcs["down"],
                            "max_fc": None,
                            "evidence": pcs["evidence"],
                            "protein_level": True,
                        }

    # Also add the global "Normal/healthy" as a reference baseline (score = 0)
    # and make sure all conditions have an entry
    all_conditions = set()
    for ac, pdata in tdata["proteins"].items():
        for site in pdata["sites"]:
            all_conditions.update(site["scores"].keys())
            all_conditions.update(site.get("conditions", {}).keys())

    tdata["conditions"] = sorted(all_conditions)
    tdata["scoring_method"] = {
        "description": "Evidence-weighted literature score per site per condition",
        "up_abstract": "+1 per PubMed abstract reporting increase",
        "down_abstract": "-1 per PubMed abstract reporting decrease",
        "up_curated": "+2 per curated UniProtKB annotation reporting increase",
        "down_curated": "-2 per curated UniProtKB annotation reporting decrease",
        "fold_change": "If abstract mentions explicit fold-change, that value is used instead",
        "interpretation": "Positive score = evidence of PTM increase; negative = decrease; magnitude = evidence strength",
    }

    # Stats
    sites_with_scores = 0
    total_scored_conditions = 0
    for ac, pdata in tdata["proteins"].items():
        for site in pdata["sites"]:
            if site["scores"]:
                sites_with_scores += 1
                total_scored_conditions += len(site["scores"])

    summary = {
        "sites_with_scores": sites_with_scores,
        "total_scored_site_conditions": total_scored_conditions,
        "unique_conditions": len(all_conditions),
        "site_level_hits": len(site_cond_scores),
        "protein_level_hits": len(protein_cond_scores),
    }

    outpath = "/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/tracks.json"
    with open(outpath, "w") as f:
        json.dump(tdata, f, separators=(",", ":"))

    import os
    print(f"Done. {os.path.getsize(outpath)} bytes", file=sys.stderr)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
