"""Process reaction records, construct a graph, and rank tandem pathways."""
import argparse
import contextlib
import csv
import io
import json
import logging
from collections import Counter
from pathlib import Path
from extraction.extract_reactions import validate
from data_processing.clean_reactions import read, write, clean_samples, metric_semantics_audit
from data_processing.clean_reactions import structure, values, abnormal
from data_processing.build_graph import build_graph
from path_recommendation.recommendation import rank
from path_recommendation.scoring import single_breakdown

ROOT = Path(__file__).resolve().parent

def run(root, output, temperature_shortlist=None):
    root, output = Path(root), Path(output)
    for logger in (values.LOG, structure.log, abnormal.LOG):
        logger.setLevel(logging.ERROR)
    with contextlib.redirect_stdout(io.StringIO()):
        raw, structured, cleaned, audit, removed = clean_samples(root)
    write(output / "structured_reactions.json", structured)
    write(output / "cleaned_reactions.json", cleaned)
    write(output / "cleaning_audit.json", audit)
    metric_audit = metric_semantics_audit(structured, cleaned)
    write(output / "performance_review.json", metric_audit)
    write(output / "removed_reactions.json", removed)
    publication_graph = read(root / 'data/graph_metadata.json')
    graph = build_graph(cleaned, read(root / "data/compound_names.json"), publication_graph)
    for name, rows in graph.items():
        write(output / "graph" / (name + ".json"), rows)
    single = []
    for art in cleaned.values():
        for rr in art["reactionDetails"]:
            explained = single_breakdown(art["articleInformation"], rr)
            single.append({"source_reaction_ID": str(rr["ID"]), "doi": art["articleInformation"]["doi"],
                "modality": rr["conditions"]["reactionType"], "score": explained["total"],
                "breakdown": explained})
    write(output / "reaction_score_breakdown.json", single)
    replay_paths, replay_stats = rank(graph["reaction_nodes"], graph["compoundToReaction_relation"], temperature_shortlist)
    write(output / "path_rankings.json", replay_paths)
    summary = {"source_graph_reactions": 471436,
        "sample_articles": len(raw), "sample_reactions_by_type": dict(Counter(r["reactionType"] for r in graph["reaction_nodes"])),
        "graph_counts": {k: len(v) for k, v in graph.items()}, "removed_reactions": len(removed),
        "performance_values_for_review": len(metric_audit),
        "ranking": replay_stats}
    write(output / "run_summary.json", summary)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "top_paths.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rank", "photo_ID", "thermal_ID", "photo_doi", "thermal_doi", "A", "B", "C", "total"])
        writer.writeheader()
        for path in replay_paths[:20]:
            writer.writerow({**{k: path[k] for k in ["rank", "photo_ID", "thermal_ID", "photo_doi", "thermal_doi", "total"]},
                             **{k: round(path["breakdown"][k], 6) for k in ["A", "B", "C"]}})
    return summary

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs', help='Results directory')
    parser.add_argument('--fulltexts', type=Path, default=ROOT / 'data/fulltexts',
                        help='Directory containing the 40 article JSON files')
    parser.add_argument('--temperature-shortlist', type=int, default=None,
                        help='Keep this many pairs with the smallest temperature differences before scoring')
    args = parser.parse_args()
    if args.temperature_shortlist is not None and args.temperature_shortlist < 1:
        parser.error('--temperature-shortlist must be positive')
    article_checks = validate(ROOT, args.fulltexts, args.output)
    summary = run(ROOT, args.output, args.temperature_shortlist)
    summary['fulltext_inputs_checked'] = article_checks['articles']
    print(json.dumps(summary, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
