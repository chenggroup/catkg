#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'Remove reactions that fail the supplied checks from flattened reactionDetails records.'

from __future__ import annotations
import json, sys, logging
from pathlib import Path
from typing import Dict, List

# Configure logging to stdout.
LOG = logging.getLogger("abnorm")
LOG.setLevel(logging.INFO)
LOG.addHandler(logging.StreamHandler(sys.stdout))


# Value helpers.
def is_unknown(val) -> bool:
    'Test for unknown strings, ignoring case and surrounding whitespace.'
    return isinstance(val, str) and val.strip().lower() == "unknown"


def norm_lower(val) -> str:
    'Strip and lowercase strings; return an empty string for other types.'
    return val.strip().lower() if isinstance(val, str) else ""


# Check one reaction.
def is_abnormal_reaction(reaction: Dict, article_id: str, idx: int,
                         rm_log: List[Dict]) -> bool:
    reactants = reaction.get("reactants", {}).get("reactantList", [])
    products  = reaction.get("products",  [])

    r_names = [norm_lower(r.get("name")) for r in reactants if r.get("name")]
    p_names = [norm_lower(p.get("name")) for p in products  if p.get("name")]

    def mark(reason: str) -> bool:
        LOG.info("Removed  %s  #%d  – %s", article_id, idx, reason)
        rm_log.append({"article_id": article_id,
                       "reaction_index": idx,
                       "reason": reason})
        return True

    # Check required fields.
    if not r_names or all(x in {"unknown", "none"} for x in r_names):
        return mark("Missing or unknown reactants")
    if not p_names or all(x in {"unknown", "none"} for x in p_names):
        return mark("Missing or unknown products")

    # Handle missing or non-string catalyst values.
    cat_name = reaction.get("catalyst", {}).get("catalystName")
    if is_unknown(cat_name):
        return mark("Catalyst name is unknown")

    # Apply reaction checks.
    ch_comp = {"ch4", "c2h2", "c2h4", "c2h6",
               "c3h8", "methanol", "ch3oh"}
    h2_src  = {"co", "h2o", "ch4", "methanol", "co2"}

    if r_names == ["co2"] and "co" not in p_names:
        return mark("CO2 -> no CO")
    if r_names == ["co"] and "co2" not in p_names:
        return mark("CO -> no CO2")
    if r_names == ["o2"]:
        return mark("Only reactant O2")
    if p_names == ["h2"] and not any(r in h2_src for r in r_names):
        return mark("Only H2 product w/o source")
    if r_names == ["h2"]:
        return mark("Only reactant H2")
    if p_names in (["o2"], ["h2o"]):
        return mark(f"Only product {p_names[0].upper()}")
    if r_names == ["h2o"]:
        return mark("Only reactant H2O")
    if r_names == ["h2o2"] and set(p_names) - {"o2", "h2o"}:
        return mark("H2O2 produces other compounds")
    if r_names == ["n2"] or p_names == ["n2"]:
        return mark("Contains only N2")
    if "f2" in r_names or "f2" in p_names:
        return mark("Contains F2")
    if r_names == ["si"]:
        return mark("Only Si reactant")

    ch_oxid_prod = {
        "co", "co2", "h2o", "c2h4o", "c2h4o2", "c2h6o",
        "ch4o", "ch2o", "acetaldehyde", "acetic acid",
        "ethanol", "methanol", "formaldehyde"
    }
    if len(r_names) == 1 and r_names[0] in ch_comp \
            and any(p in ch_oxid_prod for p in p_names):
        return mark("CH compound unexpected oxidation")

    return False


# Process the reaction collection.
def process(in_path: Path, out_path: Path, rm_log_path: Path) -> None:
    data = json.loads(in_path.read_text(encoding="utf-8"))
    rm_log: List[Dict] = []

    total_before, total_after = 0, 0

    for art_id, art in data.items():
        details = art.get("reactionDetails", [])
        total_before += len(details)

        new_details = []
        for idx, det in enumerate(details, start=1):
            if isinstance(det, dict) and not is_abnormal_reaction(det, art_id, idx, rm_log):
                new_details.append(det)

        art["reactionDetails"] = new_details
        total_after += len(new_details)

    # Write JSON and log output.
    out_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    rm_log_path.write_text(json.dumps(rm_log, ensure_ascii=False, indent=2), encoding="utf-8")

    # Report counts.
    LOG.info("Cleaned JSON written to: %s", out_path)
    LOG.info("Removed reactions log   : %s", rm_log_path)
    LOG.info("Total reactions processed: %d", total_before)
    LOG.info("Total abnormal removed   : %d", len(rm_log))
    LOG.info("Remaining valid reactions: %d", total_after)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: filter_reactions.py <input.json> <output.json> <removed_log.json>")
        sys.exit(1)
    in_f, out_f, log_f = map(Path, sys.argv[1:4])
    process(in_f, out_f, log_f)
