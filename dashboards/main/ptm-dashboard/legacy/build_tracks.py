#!/usr/bin/env python3
"""
Build a per-protein track data file for the IGV-style PTM viewer.
For each protein: sequence length, all PTM sites with positions, and
for each site: which conditions show up/down, with evidence.
"""
import json, sys, re, requests, time

def main():
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/glygen_data.json") as f:
        gdata = json.load(f)
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/quantitative.json") as f:
        qdata = json.load(f)
    with open("/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/site_evidence.json") as f:
        sdata = json.load(f)

    # Index quantitative hits by protein
    quant_by_protein = {}
    for h in qdata["hits"]:
        ac = h["protein_ac"]
        if ac not in quant_by_protein:
            quant_by_protein[ac] = []
        quant_by_protein[ac].append(h)

    # Index site evidence by protein
    sites_by_protein = {}
    for s in sdata["sites"]:
        ac = s["uniprot_ac"]
        if ac not in sites_by_protein:
            sites_by_protein[ac] = []
        sites_by_protein[ac].append(s)

    # Fetch sequences for all proteins
    proteins = {}
    for p in gdata["proteins"]:
        ac = p["uniprot_ac"]
        print(f"Fetching sequence for {ac}...", file=sys.stderr)
        try:
            r = requests.get(f"https://api.glygen.org/protein/detail/{ac}/", timeout=30)
            if r.status_code == 200:
                raw = r.json()
                seq = raw.get("sequence", {}).get("sequence", "")
                length = raw.get("sequence", {}).get("length", len(seq))
            else:
                seq = ""
                length = p.get("length", 0)
        except:
            seq = ""
            length = p.get("length", 0)

        # Collect all known PTM sites for this protein
        all_sites = {}  # pos -> site info

        # Glycosylation sites
        for g in p.get("glycosylation_sites", []):
            pos = g.get("start_pos")
            if not pos or not isinstance(pos, int):
                continue
            key = pos
            if key not in all_sites:
                all_sites[key] = {
                    "pos": pos,
                    "residue": g.get("residue", ""),
                    "site_lbl": g.get("site_lbl", f"pos{pos}"),
                    "ptm_types": set(),
                    "glytoucan": set(),
                    "conditions": {},  # condition -> {direction, evidence}
                    "tissues_context": set(),
                }
            all_sites[key]["ptm_types"].add(f"Glycosylation ({g.get('type', '')})")
            if g.get("glytoucan_ac"):
                all_sites[key]["glytoucan"].add(g["glytoucan_ac"])

        # Phosphorylation sites
        for ph in p.get("phosphorylation_sites", []):
            pos = ph.get("start_pos")
            if not pos or not isinstance(pos, int):
                continue
            key = pos
            if key not in all_sites:
                all_sites[key] = {
                    "pos": pos,
                    "residue": ph.get("residue", ""),
                    "site_lbl": ph.get("site_lbl", f"pos{pos}"),
                    "ptm_types": set(),
                    "glytoucan": set(),
                    "conditions": {},
                    "tissues_context": set(),
                }
            comment = ph.get("comment", "")
            all_sites[key]["ptm_types"].add(f"Phosphorylation")
            if comment:
                all_sites[key]["ptm_types"].add(f"Phospho: {comment[:40]}")

        # Map site evidence (tissue/condition context from abstracts)
        for se in sites_by_protein.get(ac, []):
            pos = se.get("start_pos")
            if not pos or not isinstance(pos, int):
                continue
            if pos in all_sites:
                for c in se.get("conditions", []):
                    all_sites[pos]["conditions"].setdefault(c, {"direction": "observed", "evidence": [], "count": 0})
                    all_sites[pos]["conditions"][c]["count"] += 1
                for t in se.get("tissues", []):
                    all_sites[pos]["tissues_context"].add(t)

        # Map quantitative hits (direction up/down) to specific sites
        for qh in quant_by_protein.get(ac, []):
            # Try to match sites from the quantitative hit to positions
            for site_str in qh.get("sites", []):
                # Parse site string like "Thr743", "S396", "Ser-198"
                m = re.search(r'(\d+)', site_str)
                if not m:
                    continue
                pos = int(m.group(1))
                if pos in all_sites:
                    for cond in qh.get("diseases", []):
                        entry = all_sites[pos]["conditions"].setdefault(cond, {"direction": "observed", "evidence": [], "count": 0})
                        entry["direction"] = qh["direction"]
                        entry["evidence"].append({
                            "sentence": qh["sentence"][:150],
                            "pmid": qh.get("pmid", ""),
                            "confidence": qh.get("confidence", "medium"),
                            "quantitative": qh.get("quantitative"),
                        })
                        entry["count"] += 1
                    for t in qh.get("tissues", []):
                        all_sites[pos]["tissues_context"].add(t)

            # If no specific sites, mark as protein-level for all conditions
            if not qh.get("sites"):
                for cond in qh.get("diseases", []):
                    # Store as protein-level hit (pos=0)
                    if 0 not in all_sites:
                        all_sites[0] = {
                            "pos": 0,
                            "residue": "",
                            "site_lbl": "protein-level",
                            "ptm_types": set(),
                            "glytoucan": set(),
                            "conditions": {},
                            "tissues_context": set(),
                        }
                    for pt in qh.get("ptm_types", []):
                        all_sites[0]["ptm_types"].add(pt)
                    entry = all_sites[0]["conditions"].setdefault(cond, {"direction": "observed", "evidence": [], "count": 0})
                    entry["direction"] = qh["direction"]
                    entry["count"] += 1

        # Convert sets to lists for JSON
        site_list = []
        for pos in sorted(all_sites.keys()):
            s = all_sites[pos]
            if pos == 0:
                continue  # Skip protein-level for track view
            site_list.append({
                "pos": s["pos"],
                "residue": s["residue"],
                "label": s["site_lbl"],
                "ptm_types": sorted(s["ptm_types"]),
                "glytoucan": sorted(s["glytoucan"])[:3],
                "tissues": sorted(s["tissues_context"]),
                "conditions": {
                    c: {
                        "direction": v["direction"],
                        "count": v["count"],
                        "evidence": v["evidence"][:2],
                    }
                    for c, v in s["conditions"].items()
                },
            })

        # Collect all conditions across all sites for this protein
        all_conditions = set()
        for s in site_list:
            all_conditions.update(s["conditions"].keys())

        # Diseases at protein level
        disease_names = [d["name"] for d in p.get("diseases", [])]

        proteins[ac] = {
            "ac": ac,
            "name": p["display_name"],
            "gene": p.get("gene_name", ""),
            "protein_name": p.get("protein_name", ""),
            "length": length,
            "seq_snippet": seq[:50] + "..." if len(seq) > 50 else seq,
            "sites": site_list,
            "total_sites": len(site_list),
            "conditions": sorted(all_conditions),
            "diseases": disease_names,
        }

        time.sleep(0.3)

    # Global condition list
    all_conds = set()
    for p in proteins.values():
        all_conds.update(p["conditions"])

    output = {
        "conditions": sorted(all_conds),
        "proteins": proteins,
    }

    outpath = "/app/workspace/5bd6a31a-1061-4020-a3f6-35681df676ac/ptm-dashboard/tracks.json"
    with open(outpath, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    import os
    print(f"Done. {os.path.getsize(outpath)} bytes, {len(proteins)} proteins", file=sys.stderr)
    print(json.dumps({"proteins": len(proteins), "conditions": len(all_conds), "total_sites": sum(p["total_sites"] for p in proteins.values())}, indent=2))


if __name__ == "__main__":
    main()
