#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'Apply conversion, yield, and selectivity bounds checks to flattened reactionDetails records.'

import json, re, sys
from pathlib import Path


# Percentage helpers.

def parse_pct(s: str):
    """'95%' -> 95.0 | 'unknown'/None -> None"""
    if s is None:
        return None
    s = str(s).strip()
    if s.lower() == "unknown":
        return None
    s = s.replace('\u200b', '')            # Remove zero-width characters.
    s = re.sub(r'%+', '', s).strip()
    try:
        return float(s)
    except ValueError:
        return None


def fmt_pct(val: float) -> str:
    'Format a percentage with at most two decimal places and no trailing zeros.'
    return f"{val:.2f}".rstrip('0').rstrip('.') + '%'


# Check one flattened reaction.

def filter_rr(rr: dict, rid: str):

    react_list = rr.get("reactants", {}).get("reactantList", [])
    max_conv = None
    for r in react_list:
        conv = parse_pct(r.get("conversion"))
        if conv is not None:
            if conv > 100:
                print(f"[RID {rid}] Reactant '{r.get('name','?')}' conversion {conv}% >100 -> unknown")
                r["conversion"] = "unknown"
            else:
                max_conv = max(max_conv, conv) if max_conv is not None else conv


    products = rr.get("products", [])
    tot_yld = tot_sel = 0.0
    prod_info = []  # (dict, yld, sel)
    for p in products:
        yld = parse_pct(p.get("yield"))
        sel = parse_pct(p.get("selectivity"))
        if yld is not None and yld > 100:
            print(f"[RID {rid}] Prod '{p.get('name','?')}' yield {yld}% >100 -> unknown")
            p["yield"] = "unknown"
            yld = None
        if sel is not None and sel > 100:
            print(f"[RID {rid}] Prod '{p.get('name','?')}' selectivity {sel}% >100 -> unknown")
            p["selectivity"] = "unknown"
            sel = None
        if yld is not None:
            tot_yld += yld
        if sel is not None:
            tot_sel += sel
        prod_info.append((p, yld, sel))

    # Check totals across products.
    if tot_yld > 100:
        for p, _, _ in prod_info:
            p["yield"] = "unknown"
    if tot_sel > 100:
        for p, _, _ in prod_info:
            p["selectivity"] = "unknown"

    # Apply the implemented upper-bound correction.
    for p, yld, sel in prod_info:
        if yld is None:
            continue
        if sel is not None and max_conv is not None:
            ymax = sel * max_conv / 100
            if yld * 0.95 > ymax:
                print(f"[RID {rid}] Prod '{p.get('name','?')}' yield {yld}% > theo {ymax:.2f}% -> adjusted")
                p["yield"] = fmt_pct(ymax)
        elif max_conv is not None and yld * 0.95 > max_conv:
            print(f"[RID {rid}] Prod '{p.get('name','?')}' yield {yld}% > conv {max_conv}% -> unknown")
            p["yield"] = "unknown"
        elif sel is not None and yld * 0.95 > sel:
            print(f"[RID {rid}] Prod '{p.get('name','?')}' yield {yld}% > sel {sel}% -> unknown")
            p["yield"] = "unknown"


# Traverse the JSON input.

def walk(node, prefix=""):
    if isinstance(node, dict):
        # reactionDetails contains flattened reaction dictionaries.
        if "reactionDetails" in node and isinstance(node["reactionDetails"], list):
            for rr in node["reactionDetails"]:
                if isinstance(rr, dict) and "ID" in rr:
                    rid = f"{prefix}{rr['ID']}"
                    filter_rr(rr, rid)
        for k, v in node.items():
            new_prefix = f"{prefix}{k}." if prefix or k else prefix
            walk(v, new_prefix)
    elif isinstance(node, list):
        for item in node:
            walk(item, prefix)


# Command-line entry point.

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-") and "=" not in a]
    infile = args[0] if args else r"*_test_cleaned.json"
    in_path = Path(infile)
    data = json.loads(in_path.read_text(encoding="utf-8"))
    walk(data)
    out_path = in_path.with_name(in_path.stem + "_filter.json")
    out_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print("[OK] wrote", out_path)

if __name__ == "__main__":
    main()
