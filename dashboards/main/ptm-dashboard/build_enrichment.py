#!/usr/bin/env python3
"""
Fetch GO annotations and pathway data from GlyGen for all proteins,
then compute enrichment of GO terms and pathways for PTM-affected
proteins in each disease condition (Fisher's exact test).
"""
import json, sys, time, math, requests
from collections import defaultdict

def fetch_protein_annotations(ac):
    """Fetch GO and pathway annotations from GlyGen."""
    try:
        r = requests.get(f"https://api.glygen.org/protein/detail/{ac}/", timeout=60)
        if r.status_code != 200:
            return None
        raw = r.json()

        go_terms = []
        go_data = raw.get("go_annotation", {})
        for cat in go_data.get("categories", []):
            cat_name = cat.get("name", "")
            for ann in cat.get("go_terms", []):
                go_terms.append({
                    "id": ann.get("id", ""),
                    "name": ann.get("name", ""),
                    "category": cat_name,
                    "evidence": [e.get("code", "") for e in ann.get("evidence", [])],
                })

        pathways = []
        for pw in raw.get("pathway", []):
            pathways.append({
                "id": pw.get("id", ""),
                "name": pw.get("name", ""),
                "database": pw.get("resource", ""),
            })

        return {"go_terms": go_terms, "pathways": pathways}
    except Exception as e:
        print(f"  Error fetching {ac}: {e}", file=sys.stderr)
        return None


def fisher_exact_pvalue(a, b, c, d):
    """
    One-tailed Fisher's exact test p-value (enrichment).
    Contingency table:
        | In set | Not in set |
    Has term  |   a   |     b     |
    No term   |   c   |     d     |
    """
    # Use logarithmic calculation for factorials to avoid overflow
    def log_factorial(n):
        return sum(math.log(i) for i in range(1, n + 1)) if n > 0 else 0

    n = a + b + c + d
    # P-value = sum of probabilities of tables as or more extreme than observed
    # For one-tailed (enrichment): fix margins, sum over a' >= a
    row1 = a + b
    col1 = a + c
    row2 = c + d
    col2 = b + d

    log_denom = log_factorial(n)
    log_margins = (log_factorial(row1) + log_factorial(row2) +
                   log_factorial(col1) + log_factorial(col2) - log_denom)

    p_value = 0.0
    max_a = min(row1, col1)
    for ai in range(a, max_a + 1):
        bi = row1 - ai
        ci = col1 - ai
        di = row2 - ci
        if bi < 0 or ci < 0 or di < 0:
            continue
        log_p = log_margins - (log_factorial(ai) + log_factorial(bi) +
                               log_factorial(ci) + log_factorial(di))
        p_value += math.exp(log_p)

    return min(p_value, 1.0)


def main():
    # Load existing data
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/glygen_data.json") as f:
        gdata = json.load(f)
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/quantitative.json") as f:
        qdata = json.load(f)

    proteins = {p["uniprot_ac"]: p for p in gdata["proteins"]}
    all_acs = list(proteins.keys())

    # Fetch GO and pathway annotations
    print("Fetching GO/pathway annotations...", file=sys.stderr)
    annotations = {}
    for i, ac in enumerate(all_acs):
        print(f"  [{i+1}/{len(all_acs)}] {ac}", file=sys.stderr)
        ann = fetch_protein_annotations(ac)
        if ann:
            annotations[ac] = ann
        time.sleep(0.4)

    # Build term-to-protein mappings
    # GO terms
    go_to_proteins = defaultdict(set)  # go_id -> set of protein ACs
    protein_go = defaultdict(set)       # ac -> set of go_ids
    go_info = {}                        # go_id -> {name, category}

    for ac, ann in annotations.items():
        for gt in ann["go_terms"]:
            gid = gt["id"]
            go_to_proteins[gid].add(ac)
            protein_go[ac].add(gid)
            if gid not in go_info:
                go_info[gid] = {"name": gt["name"], "category": gt["category"]}

    # Pathways
    pw_to_proteins = defaultdict(set)
    protein_pw = defaultdict(set)
    pw_info = {}

    for ac, ann in annotations.items():
        for pw in ann["pathways"]:
            pid = pw["id"]
            pw_to_proteins[pid].add(ac)
            protein_pw[ac].add(pid)
            if pid not in pw_info:
                pw_info[pid] = {"name": pw["name"], "database": pw["database"]}

    print(f"GO terms: {len(go_info)}, Pathways: {len(pw_info)}", file=sys.stderr)

    # Build condition -> affected proteins map from quantitative data
    cond_proteins = defaultdict(set)  # condition -> set of protein ACs with PTM changes
    cond_protein_direction = defaultdict(lambda: defaultdict(list))  # condition -> ac -> [directions]

    for hit in qdata["hits"]:
        ac = hit["protein_ac"]
        for cond in hit.get("diseases", []):
            if cond in ("Control", "Normal/healthy"):
                continue
            cond_proteins[cond].add(ac)
            cond_protein_direction[cond][ac].append(hit["direction"])

    print(f"Conditions with affected proteins: {len(cond_proteins)}", file=sys.stderr)

    # Compute enrichment for each condition
    N = len(all_acs)  # background size

    enrichment_results = {}

    for cond, affected_set in cond_proteins.items():
        if len(affected_set) < 2:
            continue  # Need at least 2 proteins for meaningful enrichment

        k = len(affected_set)  # size of test set

        # GO enrichment
        go_results = []
        for gid, gproteins in go_to_proteins.items():
            m = len(gproteins)  # proteins with this term in background
            x = len(affected_set & gproteins)  # overlap
            if x == 0:
                continue

            # Fisher's exact test
            a = x           # affected AND has term
            b = m - x       # not affected AND has term
            c = k - x       # affected AND no term
            d = N - m - c   # not affected AND no term
            if d < 0:
                d = 0

            pval = fisher_exact_pvalue(a, b, c, d)
            fold = (x / k) / (m / N) if (m / N) > 0 else 0

            go_results.append({
                "id": gid,
                "name": go_info[gid]["name"],
                "category": go_info[gid]["category"],
                "overlap": x,
                "term_size": m,
                "set_size": k,
                "bg_size": N,
                "fold_enrichment": round(fold, 2),
                "p_value": round(pval, 6),
                "proteins": sorted(affected_set & gproteins),
            })

        # Pathway enrichment
        pw_results = []
        for pid, pproteins in pw_to_proteins.items():
            m = len(pproteins)
            x = len(affected_set & pproteins)
            if x == 0:
                continue

            a = x
            b = m - x
            c = k - x
            d = N - m - c
            if d < 0:
                d = 0

            pval = fisher_exact_pvalue(a, b, c, d)
            fold = (x / k) / (m / N) if (m / N) > 0 else 0

            pw_results.append({
                "id": pid,
                "name": pw_info[pid]["name"],
                "database": pw_info[pid]["database"],
                "overlap": x,
                "term_size": m,
                "set_size": k,
                "bg_size": N,
                "fold_enrichment": round(fold, 2),
                "p_value": round(pval, 6),
                "proteins": sorted(affected_set & pproteins),
            })

        # Sort by p-value
        go_results.sort(key=lambda x: x["p_value"])
        pw_results.sort(key=lambda x: x["p_value"])

        # Protein details for this condition
        prot_details = []
        for ac in sorted(affected_set):
            pname = proteins[ac]["display_name"] if ac in proteins else ac
            dirs = cond_protein_direction[cond][ac]
            net_dir = "up" if dirs.count("up") > dirs.count("down") else "down" if dirs.count("down") > dirs.count("up") else "mixed"
            prot_details.append({"ac": ac, "name": pname, "direction": net_dir, "hits": len(dirs)})

        enrichment_results[cond] = {
            "affected_proteins": len(affected_set),
            "protein_details": prot_details,
            "go_enrichment": go_results[:50],  # top 50
            "pathway_enrichment": pw_results[:30],  # top 30
        }

    # Also save the full GO/pathway annotations per protein for reference
    protein_annotations = {}
    for ac in all_acs:
        ann = annotations.get(ac, {"go_terms": [], "pathways": []})
        protein_annotations[ac] = {
            "name": proteins[ac]["display_name"] if ac in proteins else ac,
            "go_count": len(ann["go_terms"]),
            "pathway_count": len(ann["pathways"]),
            "go_categories": list(set(g["category"] for g in ann["go_terms"])),
            "pathways": [{"id": p["id"], "name": p["name"]} for p in ann["pathways"]],
        }

    output = {
        "summary": {
            "total_go_terms": len(go_info),
            "total_pathways": len(pw_info),
            "conditions_analyzed": len(enrichment_results),
            "background_size": N,
        },
        "enrichment": enrichment_results,
        "protein_annotations": protein_annotations,
    }

    outpath = "/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/enrichment.json"
    with open(outpath, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    import os
    print(f"\nDone. {os.path.getsize(outpath)} bytes", file=sys.stderr)
    print(json.dumps(output["summary"], indent=2))


if __name__ == "__main__":
    main()
