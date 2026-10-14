"""Print the scoring calculation for the highest-ranked sample path.

Run: python path_recommendation/scoring_and_ranking_example.py
Use --source neo4j to query an existing graph and recalculate its sample paths.
For the supplied sample, the highest-scoring path is reaction 284292 -> 311263.
This example explains that path; it does not imply experimental validation.
"""
import argparse
import getpass
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from data_processing.clean_reactions import write
from path_recommendation.recommendation import (
    load_ranked_paths, retrieve_and_rank, sample_reaction_ids, compound_names, formatted_rankings,
)
from path_recommendation.scoring import SIMPLE, as_list, path_breakdown, single_breakdown

def show_reaction_score(step, node, info, reaction):
    """Print A1-A8 using the existing single-reaction scoring functions."""
    detail = single_breakdown(info, reaction)
    terms = detail["terms"]
    conditions = reaction["conditions"]
    print(f"\nStep {step}: reaction {node['ID']} ({node['reactionType']})")
    print(f"  DOI: {info['doi']}")
    print(f"  Catalyst: {reaction['catalyst']['catalystName']}")

    impact = float(info.get("impactEffect", 0) or 0)
    print(f"  A1 Journal impact: round({impact:g} / 50.5 * 200, 2) = {terms['A1_journal_impact']:.2f}")
    print(f"  A2 Publication date: {info.get('pubDate')} -> {terms['A2_publication_date']:.2f}")
    print("     The date score increases linearly from 2000-01-01 to 2024-07-30.")
    print(f"  A3 Citations: {info.get('citations', 0)} -> {terms['A3_citations']:.2f}")
    print("     Counts up to 100 contribute their value; counts above 100 contribute 200.")
    fields = detail["known_catalyst_fields"]
    print(f"  A4 Catalyst fields: {len(fields)} * 20 = {terms['A4_catalyst_completeness']:.2f}")
    print(f"     Available fields: {', '.join(fields) or 'none'}")
    print(f"  A5 Temperature available: {conditions.get('temperature')} -> {terms['A5_temperature_available']:.2f}")
    print(f"  A6 Gas names and pressures available: {terms['A6_all_gas_pressures_available']:.2f}")
    for gas in conditions.get("gases", []):
        print(f"     {gas.get('name')}: {gas.get('pressure', 'unknown')}")
    print(f"  A7 Reaction time available: {conditions.get('reactionTime')} -> {terms['A7_reaction_time_available']:.2f}")

    conversions = reaction["reactants"]["reactantList"]
    print("  A8 Product performance:")
    print("     Conversion inputs: " + ", ".join(
        f"{item['name']}={item.get('conversion', 'unknown')}" for item in conversions
    ))
    print("     x = maximum reactant conversion; s = selectivity; y = yield (fractions).")
    print("     Missing percentages use 0 in the calculation; all recorded products contribute.")
    for product in detail["products"]:
        x = product["max_conversion_fraction"]
        s = product["selectivity_fraction"]
        y = product["yield_fraction"]
        print(f"     {product['name']}: yield={product['yield_cleaned']}, selectivity={product['selectivity_cleaned']}")
        print(f"       400 * max({s:g} * {x:g}, {y:g}) + 100 * (2 * {s:g} + {x:g} + {y:g})")
        print(f"       = {product['contribution']:.6f}")
    print(f"     A8 = sum of product contributions = {terms['A8_performance']:.6f}")
    print("  A_step = round(A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8, 2)")
    print("         = round(" + " + ".join(f"{value:.6f}" for value in terms.values()) + ", 2)")
    print(f"         = {detail['total']:.2f}")

    # Confirm that the explanation agrees with scores assigned during cleaning.
    if any(not math.isclose(float(product["score"]), detail["total"], abs_tol=1e-8)
           for product in reaction["products"]):
        raise AssertionError("The displayed score differs from the assigned reaction score.")
    return detail["total"]

def show_pair_score(forward_steps, detail):
    """Print the input comparisons and B1-B8 contributions for adjacent steps."""
    first, second = forward_steps
    terms = detail["terms"]
    print("\nB: compatibility of the two adjacent reactions")
    for code, label, field in [
        ("B1", "Heterogeneous status", "isHeterogeneous"),
        ("B2", "Catalyst name", "catalystName"),
        ("B3", "Catalyst support", "catalystSupport"),
    ]:
        print(f"  {code} {label}: {first.get(field, 'unknown')} / {second.get(field, 'unknown')} -> {terms[code]:.2f}")
    promoters = [sorted(item.get("promoterName", "") for item in as_list(step.get("promoters", [])))
                 for step in forward_steps]
    print(f"  B4 Promoter sets: {promoters[0]} / {promoters[1]} -> {terms['B4']:.2f}")
    difference = detail["temperature_difference_K"]
    if difference is None:
        print(f"  B5 Temperature: missing numeric temperatures -> {terms['B5']:.2f}")
    else:
        print(f"  B5 Temperature: {first['temperature']} / {second['temperature']}")
        print(f"     30000 / ({difference:g} + 30) = {terms['B5']:.2f}")
    gases = [as_list(step.get("gases", [])) for step in forward_steps]
    gas_names = [sorted(gas.get("name", "") for gas in group) for group in gases]
    print(f"  B6 Shared gas names: {gas_names[0]} / {gas_names[1]} -> {terms['B6']:.2f}")
    print(f"  B7 Shared-gas pressure agreement: {terms['B7']:.2f}")
    if terms["B6"] == 0:
        print("     No shared gas name, so no pressure comparison contributes.")
    print("     The largest usable shared-gas score is 200 * (1 - abs(p1-p2)/(p1+p2)).")
    solvents = [as_list(step.get("solventList", [])) for step in forward_steps]
    print(f"  B8 Shared solvent names: {solvents[0]} / {solvents[1]} -> {terms['B8']:.2f}")
    print("  B = " + " + ".join(f"{terms[f'B{number}']:.2f}" for number in range(1, 9))
          + f" = {detail['total']:.2f}")

def show_ranked_paths(ranked, records, source="the supplied sample"):
    """Print the calculation for the highest-ranked path from either data source."""
    selected, score = ranked[0]
    detail = path_breakdown(selected)
    # Queries store steps in backward order; display them in reaction order.
    forward_steps = list(reversed(selected["reaction"]))
    identifiers = " -> ".join(str(step["ID"]) for step in forward_steps)

    print("Scoring and ranking example")
    print(f"Ranked candidate paths: {len(ranked)}")
    print(f"Selected reaction IDs: {identifiers}")
    print(f"This demonstration explains only the highest-scoring path in {source}.")
    print("Reaction sequence: CO2 + H2 -> CO; CO + H2 -> methanol.")
    print("\nFirst five paths ranked by the project scoring function:")
    for position, (path, total) in enumerate(ranked[:5], 1):
        ids = " -> ".join(str(step["ID"]) for step in reversed(path["reaction"]))
        print(f"  {position:2d}. {ids}: {total:.2f}")

    step_scores = []
    for number, node in enumerate(forward_steps, 1):
        info, reaction = records[str(node["ID"])]
        step_scores.append(show_reaction_score(number, node, info, reaction))
    if not math.isclose(sum(step_scores), detail["A"], abs_tol=1e-8):
        raise AssertionError("The reaction scores differ from the path's Category A contribution.")
    print("\nA: sum of the single-reaction scores")
    print("  A = " + " + ".join(f"{value:.2f}" for value in step_scores) + f" = {detail['A']:.2f}")
    show_pair_score(forward_steps, detail["B_pairs"][0])

    print("\nC: whole-path contributions")
    feed = as_list(forward_steps[0].get("reactantList", []))
    print(f"  C1 First-step reactants: {feed}")
    print(f"     All are in the simple-compound list: {all(name in SIMPLE for name in feed)} -> {detail['C1_simple_feed']:.2f}")
    compounds = detail["C2_counted_compounds"]
    print(f"  C2 Counted compound names ({len(compounds)}): {', '.join(compounds)}")
    if len(compounds) <= 15:
        print(f"     10 * (15 - {len(compounds)}) = {detail['C2_fewer_compounds']:.2f}")
    else:
        print("     More than 15 names: C2 = 0.00")
    print("     This term counts stored names, including separately recorded aliases.")
    print(f"  C = {detail['C1_simple_feed']:.2f} + {detail['C2_fewer_compounds']:.2f} = {detail['C']:.2f}")

    print("\nTotal: round(A + B + C, 2)")
    print(f"  round({detail['A']:.2f} + {detail['B']:.2f} + {detail['C']:.2f}, 2) = {detail['total']:.2f}")
    if not math.isclose(detail["total"], score, abs_tol=1e-8):
        raise AssertionError("The explanation differs from the pathway scoring function.")
    print(f"Final score: {score:.2f}")
    print("Verified against calculate_score_between_reactions().")

def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("offline", "neo4j"), default="offline")
    parser.add_argument("--uri", default=os.getenv("NEO4J_URI", "bolt://localhost:7687"))
    parser.add_argument("--database", default=os.getenv("NEO4J_DATABASE", "neo4j"))
    parser.add_argument("--username", default=os.getenv("NEO4J_USERNAME", "neo4j"))
    parser.add_argument("--no-auth", action="store_true", help="Use only if database authentication is disabled.")
    parser.add_argument("--output", type=Path, help="Save all rankings and a query summary to this directory.")
    args = parser.parse_args()
    if args.source == "offline":
        ranked, records = load_ranked_paths(ROOT)
        summary = {"source": "offline", "sample_reactions": len(records), "ranked_pairs": len(ranked)}
        print("Data source: cleaned sample records (offline)")
        show_ranked_paths(ranked, records)
    else:
        try:
            from neo4j import GraphDatabase, READ_ACCESS
        except ImportError:
            parser.error("Install the Neo4j driver: python -m pip install neo4j==5.28.6")
        password = None if args.no_auth else os.getenv("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
        auth = None if args.no_auth else (args.username, password)
        with GraphDatabase.driver(args.uri, auth=auth, connection_timeout=10) as driver:
            with driver.session(database=args.database, default_access_mode=READ_ACCESS) as session:
                ranked, records, summary = session.execute_read(
                    retrieve_and_rank, sample_reaction_ids(ROOT), compound_names(ROOT))
        print("Data source: Neo4j (read-only queries)")
        print(f"Database reaction count: {summary['database_reactions']}")
        print(f"Sample reactions read from Neo4j: {summary['sample_reactions']}")
        print(f"Connected path pairs retrieved from Neo4j: {summary['retrieved_pairs']}")
        print("Single-reaction scores were recalculated from database properties.")
        print(f"Recalculated scores matching stored product scores: {summary['stored_score_matches']} reactions")
        show_ranked_paths(ranked, records, "the sample paths retrieved from Neo4j")
    output = args.output or ROOT / "outputs" / ("neo4j_scoring" if args.source == "neo4j" else "scoring_example")
    write(output / "path_rankings.json", formatted_rankings(ranked, records))
    write(output / "verification.json", summary)
    print(f"\nRanking and query summary saved to: {output.resolve()}")


if __name__ == "__main__":
    main()
