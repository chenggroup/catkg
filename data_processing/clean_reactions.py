"""Clean the saved reactions and retain their source identifiers."""
import copy
import json
import re
from pathlib import Path
from . import normalize_structure as structure, normalize_values as values
from . import check_bounds as bounds, filter_reactions as abnormal
from path_recommendation.scoring import add_scores_and_single_flag

def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")

def clean_samples(root):
    raw = read(root / "data/extracted_reactions.json")
    snapshot = read(root / "data/graph_metadata.json")
    source_flags = {r["ID"]: r["singleReaction"] for r in snapshot["reaction_nodes"]}
    substance_map = read(root / "data/compound_names.json")
    values.SUBSTANCE_MAP = {k: {"cid": v.get("cid"), "cidName": v.get("cid_name"),
                               "chebiID": v.get("chebi_id")} for k, v in substance_map.items()}
    audit = []
    # Normalize extraction fields before applying value and bounds checks.
    structured = structure.walk(copy.deepcopy(raw))
    cleaned = copy.deepcopy(structured)
    values._traverse(cleaned)
    bounds.walk(cleaned)
    removed = []
    for aid, article in cleaned.items():
        article["reactionDetails"] = [rr for i, rr in enumerate(article["reactionDetails"], 1)
            if not abnormal.is_abnormal_reaction(rr, aid, i, removed)]
    for aid, article in cleaned.items():
        original_products = structured[aid]["reactionDetails"][0]["products"]
        for reaction in article["reactionDetails"]:
            for old, product in zip(original_products, reaction["products"]):
                if old.get("yield") != product.get("yield") or old.get("selectivity") != product.get("selectivity"):
                    audit.append({"source_reaction_ID": str(reaction["ID"]), "doi": article["articleInformation"]["doi"],
                        "product": product["name"], "raw_yield": old.get("yield"), "cleaned_yield": product.get("yield"),
                        "raw_selectivity": old.get("selectivity"), "cleaned_selectivity": product.get("selectivity"),
                        "rule": "Percentage normalization and bounds checks"})
    add_scores_and_single_flag(cleaned)
    # Keep the reaction-classification flag stored in the graph.
    for article in cleaned.values():
        for reaction in article["reactionDetails"]:
            reaction["conditions"]["singleReaction"] = source_flags[str(reaction["ID"])]
    return raw, structured, cleaned, audit, removed

def metric_semantics_audit(structured, cleaned):
    """List unit-bearing performance values alongside their cleaned counterparts."""
    audit = []
    for aid, art in structured.items():
        for rr in art["reactionDetails"]:
            final = next((r for r in cleaned[aid]["reactionDetails"] if r["ID"] == rr["ID"]), None)
            if final is None:
                continue
            for product, cp in zip(rr["products"], final["products"]):
                raw_yield = str(product.get("yield", "unknown"))
                unit_signal = bool(re.search(r"(?:[μµu]?mol|\bmmol\b|\bTON\b|\bTOF\b|gcat|g-catal)", raw_yield, re.I))
                if not unit_signal:
                    continue
                audit.append({"source_reaction_ID": str(rr["ID"]), "doi": art["articleInformation"]["doi"],
                    "product": product["name"], "raw_yield": raw_yield, "cleaned_yield": cp.get("yield"),
                    "contains_percent_elsewhere": "%" in raw_yield,
                    "flag": "rate_or_amount",
                    "percent_context_check": "%" in raw_yield and cp.get("yield") != "unknown"})
    return audit
