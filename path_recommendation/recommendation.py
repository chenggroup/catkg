"""Retrieve connected sample paths offline or from Neo4j and rank them."""
import ast
import contextlib
import copy
import csv
import io
import json
import logging
import math
from pathlib import Path

from data_processing.clean_reactions import clean_samples, read
from data_processing.build_graph import build_graph, mark_main_products
from . import scoring as path_scoring
from .scoring import CATALYST_FIELDS, add_scores_and_single_flag, single_breakdown, path_breakdown
from .scoring import calculate_score_between_reactions, filter_data
from .queries import REACTION_QUERY, PATH_QUERY

ROOT = Path(__file__).resolve().parents[1]

def candidate_results(nodes, edges, require_main=True, max_products=8):
    """Retrieve CO2/H2 -> CO -> methanol through actual compound edges."""
    by_r = {}
    for edge in edges:
        by_r.setdefault(edge["r"], []).append(edge)
    groups = {"photocatalysis": [], "thermal catalysis": []}
    for node in sorted(nodes, key=lambda r: int(r["ID"])):
        mode = node.get("reactionType")
        if mode not in groups:
            continue
        relations = by_r.get(node["ID"], [])
        reactants = {e["c"] for e in relations if e["t"] == "reactant"}
        products = [e for e in relations if e["t"] == "product"]
        expected = {"Carbon Dioxide", "Hydrogen"} if mode == "photocatalysis" else {"Carbon Monoxide", "Hydrogen"}
        target = "Carbon Monoxide" if mode == "photocatalysis" else "Methanol"
        target_edge = next((e for e in products if e["c"] == target), None)
        try:
            temperature = float(node.get("temperature", "").replace(" K", ""))
        except ValueError:
            continue
        if not math.isfinite(temperature) or reactants != expected or len(products) > max_products or target_edge is None:
            continue
        if require_main and target_edge["p"].get("mainProduct") != "true":
            continue
        groups[mode].append((node, target_edge["p"].get("score", 0)))
    results = []
    # Backward retrieval stores the target-side reaction first.
    for photo, pscore in groups["photocatalysis"]:
        for thermal, tscore in groups["thermal catalysis"]:
            results.append({"score": [tscore, pscore], "reaction": [thermal, photo],
                            "path": ["Carbon Monoxide", ["Carbon Dioxide", "Hydrogen"]]})
    return path_scoring.filter_data(results), {k: len(v) for k, v in groups.items()}

def rank(nodes, edges, temperature_shortlist=None):
    results, groups = candidate_results(nodes, mark_main_products(copy.deepcopy(edges)))
    if temperature_shortlist is not None:
        results = path_scoring.calculate_temperature_difference(results, temperature_shortlist)
    ranked = path_scoring.calculate_score_between_reactions(results)
    output = formatted_rankings(ranked)
    return output, {"eligible_reactions": groups, "candidate_pairs_before_score_exclusion": len(results),
                    "ranked_pairs": len(output), "temperature_shortlist": temperature_shortlist,
                    "query_order": "target-side first: thermal, photo", "ranking": "A+B+C"}

def sample_reaction_ids(root):
    """The sample index defines the query scope, without supplying scoring values."""
    with (root / "data/samples.csv").open(encoding="utf-8", newline="") as handle:
        return [row["reaction_id"] for row in csv.DictReader(handle)]

def decode(value):
    """Decode lists stored as JSON or Python literals; never execute their content."""
    if isinstance(value, str) and value.strip().startswith(("[", "{")):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return ast.literal_eval(value)
    return value

def reaction_record(row):
    """Rebuild the scoring input from reaction, catalyst, article, and edge properties."""
    node = dict(row["node"])
    rid = str(node["ID"])
    if len(row["articles"]) != 1 or len(row["catalysts"]) != 1:
        raise ValueError(f"Reaction {rid} requires one linked article and one linked catalyst.")
    info = {key: node[key] for key in ("doi", "pubDate", "citations", "impactEffect") if key in node}
    info.update(row["articles"][0])
    cat = {key: decode(node[key]) for key in CATALYST_FIELDS + ["catalystName"] if key in node}
    cat.update({key: decode(value) for key, value in row["catalysts"][0].items()})
    conditions = {key: decode(value) for key, value in node.items()}
    conditions.setdefault("gases", [])

    def compounds(edges):
        return [{**edge["props"], "name": edge["props"].get("originName") or edge["name"],
                 "cidName": edge["name"]} for edge in edges if edge["name"] is not None]

    reaction = {"ID": rid, "reactionType": node.get("reactionType"), "catalyst": cat,
                "conditions": conditions, "reactants": {"reactantList": compounds(row["reactants"])},
                "products": compounds(row["products"])}
    stored_scores = [product.get("score") for product in reaction["products"]]
    add_scores_and_single_flag({"article": {"articleInformation": info, "reactionDetails": [reaction]}})
    # Preserve the database flag: a demonstration subset does not change article scope.
    reaction["conditions"]["singleReaction"] = node.get("singleReaction", "unknown")
    return node, info, reaction, stored_scores

def retrieve_and_rank(tx, reaction_ids, names):
    """Retrieve actual graph paths and score them in Python within a read transaction."""
    count = tx.run("MATCH (r:Reaction) RETURN count(r) AS n").single()["n"]
    rows = tx.run(REACTION_QUERY, reaction_ids=reaction_ids).data()
    nodes, records, scores = {}, {}, {}
    stored_match = 0
    for row in rows:
        node, info, reaction, stored = reaction_record(row)
        rid = str(node["ID"])
        nodes[rid], records[rid] = node, (info, reaction)
        scores[rid] = single_breakdown(info, reaction)["total"]
        try:
            matching = bool(stored) and all(value is not None and abs(float(value) - scores[rid]) < 1e-8 for value in stored)
        except (TypeError, ValueError):
            matching = False
        if matching:
            stored_match += 1
    missing = set(reaction_ids) - records.keys()
    if missing:
        raise ValueError("The database is missing sample reactions: " + ", ".join(sorted(missing)))

    pairs = tx.run(PATH_QUERY, reaction_ids=reaction_ids, **names).data()
    candidates = []
    for pair in pairs:
        photo, thermal = str(pair["photo_ID"]), str(pair["thermal_ID"])
        candidates.append({"reaction": [nodes[thermal], nodes[photo]],
                           "score": [scores[thermal], scores[photo]],
                           "path": [names["intermediate"], [names["feed"], names["hydrogen"]]]})
    ranked = calculate_score_between_reactions(filter_data(candidates))
    if not ranked:
        raise ValueError("No eligible connected paths were found for the sample reactions.")
    summary = {"database_reactions": count, "sample_reactions": len(records),
               "retrieved_pairs": len(pairs), "ranked_pairs": len(ranked),
               "stored_score_matches": stored_match}
    return ranked, records, summary

def compound_names(root):
    """Resolve query compounds using the supplied canonical-name mapping."""
    mapping = json.loads((root / "data/compound_names.json").read_text(encoding="utf-8"))
    return {key: mapping[name]["cid_name"] for key, name in
            [("feed", "CO2"), ("hydrogen", "H2"), ("intermediate", "CO"), ("target", "CH3OH")]}

def load_ranked_paths(root=ROOT):
    """Clean the sample and retrieve and rank paths with the project functions."""
    root = Path(root)
    loggers = (logging.getLogger('clean-val'), logging.getLogger('clean_refined'), logging.getLogger('abnorm'))
    previous_levels = [logger.level for logger in loggers]
    try:
        for logger in loggers:
            logger.setLevel(logging.ERROR)
        with contextlib.redirect_stdout(io.StringIO()):
            _, _, cleaned, _, _ = clean_samples(root)
    finally:
        for logger, level in zip(loggers, previous_levels):
            logger.setLevel(level)

    graph = build_graph(
        cleaned,
        read(root / "data/compound_names.json"),
        read(root / "data/graph_metadata.json"),
    )
    candidates, _ = candidate_results(
        graph["reaction_nodes"], graph["compoundToReaction_relation"]
    )
    ranked = calculate_score_between_reactions(candidates)
    if not ranked:
        raise RuntimeError("No eligible paths were found in the supplied sample.")

    records = {}
    for article in cleaned.values():
        for reaction in article["reactionDetails"]:
            records[str(reaction["ID"])] = (article["articleInformation"], reaction)
    return ranked, records

def formatted_rankings(ranked, records=None):
    """Export forward reaction IDs, source DOIs, and score components."""
    output = []
    for position, (path, total) in enumerate(ranked, 1):
        thermal, photo = path["reaction"]
        output.append({"rank": position, "photo_ID": photo["ID"], "thermal_ID": thermal["ID"],
                       "photo_doi": records[str(photo["ID"])][0]["doi"] if records else photo["doi"],
                       "thermal_doi": records[str(thermal["ID"])][0]["doi"] if records else thermal["doi"],
                       "total": total, "breakdown": path_breakdown(path)})
    return output
