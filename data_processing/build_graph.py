"""Construct reaction, compound, catalyst, and publication records."""
from __future__ import annotations
import json
import copy
from collections import Counter
from typing import Any, Dict, List


# Adapted from CatKG_project/code/clean/neo4j_article_reaction_nodes.py
def extract_reactions(data: Dict[str, Any], start_id: int) -> tuple[int, List[Dict[str, Any]]]:
    next_id = start_id
    new_nodes: List[Dict[str, Any]] = []
    for article in data.values():
        info = article.get("articleInformation", {})
        title        = info.get("title")
        doi          = info.get("doi")
        pub_date     = info.get("pubDate")
        citations    = info.get("citations")
        impact       = info.get("impactEffect")
        if isinstance(pub_date, list):   pub_date  = pub_date[0] if pub_date else None
        if isinstance(citations, list):  citations = citations[0] if citations else None

        for reaction in article.get("reactionDetails", []):
            rid = str(next_id); next_id += 1
            reaction["ID"] = rid

            node = {"ID": rid}
            if title:        node["title"]        = title
            if doi:          node["doi"]          = doi
            if pub_date:     node["pubDate"]      = pub_date
            if citations is not None: node["citations"] = citations
            if impact:       node["impactEffect"] = impact

            # Catalyst
            cat = reaction.get("catalyst", {})
            for k in ("catalystName","catalystSupport","catPreparationMethod"):
                v = cat.get(k)
                if isinstance(v, str) and v.lower() != "unknown": node[k] = v
                elif v: node[k] = v
            if cat.get("promoters"):
                node["promoters"] = json.dumps(cat["promoters"], ensure_ascii=False)

            # Reactants
            react_strs, react_names = [], []
            for r in reaction.get("reactants", {}).get("reactantList", []):
                n, conv = r.get("cidName") or r.get("name"), r.get("conversion")
                if isinstance(n, str) and n.lower() != "unknown":
                    react_names.append(n)
                    react_strs.append(f"{n}(conv.{conv})" if isinstance(conv, str) and conv.lower() != "unknown" else n)
            if react_names: node["reactantList"] = json.dumps(react_names, ensure_ascii=False)

            # Products
            prod_strs, prod_names = [], []
            for p in reaction.get("products", []):
                n  = p.get("cidName") or p.get("name")
                y  = p.get("yield");  s = p.get("selectivity")
                if isinstance(n, str) and n.lower() != "unknown":
                    prod_names.append(n)
                    label = ""
                    if isinstance(s, str) and s.lower() != "unknown": label += f"(sele.{s})"
                    if (isinstance(y, str) and y.lower() != "unknown") or isinstance(y,(int,float)): label += f"(yield.{y})"
                    prod_strs.append(f"{n}{label}")
            if prod_names: node["productList"] = json.dumps(prod_names, ensure_ascii=False)

            if react_strs and prod_strs:
                node["reactionEquation"] = " + ".join(react_strs) + " ----> " + " + ".join(prod_strs)

            # Solvents
            solv_names, solvent_ph = [], None
            for s in reaction.get("solvents", {}).get("solventList", []):
                nm, ph = s.get("cidName") or s.get("name"), s.get("solventPH")
                if isinstance(nm, str) and nm.lower() != "unknown": solv_names.append(nm)
                if isinstance(ph, str) and ph.lower() != "unknown": solvent_ph = ph
            if solv_names: node["solventList"] = json.dumps(solv_names, ensure_ascii=False)
            if solvent_ph:
                node["solventPH"] = solvent_ph
                node["solvent_pH"] = solvent_ph

            # Other conditions
            reaction_type = (
                reaction.get("reactionType")
                or reaction.get("reaction_type")
                or reaction.get("reactionClass")
            )
            if reaction_type not in (None, "", "unknown"):
                node["reactionType"] = (
                    json.dumps(reaction_type, ensure_ascii=False)
                    if isinstance(reaction_type, (list, dict))
                    else reaction_type
                )
            for k, v in reaction.get("conditions", {}).items():
                if v not in (None, "unknown", [], {}):
                    node[k] = json.dumps(v, ensure_ascii=False) if isinstance(v,(list,dict)) else v

            new_nodes.append(node)
    return next_id, new_nodes


# Adapted from CatKG_project/code/clean/neo4j_article_reaction_nodes.py
def article_node(art_id: int, info: dict) -> dict:
    get = lambda k: (info.get(k) or [""])[0] if isinstance(info.get(k), list) else info.get(k, "")
    abs_val = info.get("abstract", "")
    abstract = " ".join(abs_val) if isinstance(abs_val, list) else abs_val.strip() if isinstance(abs_val, str) else ""
    node = {
        "article_ID": str(art_id),
        "title":        get("title"),
        "journalName":  get("journalName"),
        "authors":      ", ".join(info.get("authors", [])),
        "pubDate":      get("pubDate"),
        "citations":    get("citations"),
        "doi":          get("doi"),
        "abstract":     abstract,
        "impactEffect": get("impactEffect")
    }
    return {k: v for k, v in node.items() if v}


# Adapted from CatKG_project/code/clean/relation_creat.py
def is_unknown(val: Any) -> bool:
    return isinstance(val, str) and val.strip().lower() == "unknown"


# Adapted from CatKG_project/code/clean/relation_creat.py
def _clean_name(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    name = value.strip()
    return name or None


# Adapted from CatKG_project/code/clean/relation_creat.py
def extract_compound(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    rels: List[Dict[str, Any]] = []
    for art in data.values():
        for det in art.get("reactionDetails", []):
            rid = str(det.get("ID"))
            for r in det.get("reactants", {}).get("reactantList", []):
                name = _clean_name(r.get("cidName") or r.get("name"))
                if not name:
                    continue
                origin = _clean_name(r.get("name"))
                props = {"originName": origin} if origin else {}
                c = r.get("conversion")
                if c is not None and not is_unknown(c):
                    props["conversion"] = c
                rels.append({"c": name, "r": rid, "t": "reactant", "p": props})

            for p in det.get("products", []):
                name = _clean_name(p.get("cidName") or p.get("name"))
                if not name:
                    continue
                origin = _clean_name(p.get("name"))
                props = {"originName": origin} if origin else {}
                for k in ("yield", "selectivity", "faradaicEfficiency", "score"):
                    v = p.get(k)
                    if v is not None and not is_unknown(str(v)):
                        props[k] = float(v) if k == "score" else v
                rels.append({"c": name, "r": rid, "t": "product", "p": props})

            sph = det.get("solvents", {}).get("solventPH")
            if is_unknown(sph):
                sph = None
            for s in det.get("solvents", {}).get("solventList", []):
                name = _clean_name(s.get("cidName") or s.get("name"))
                if not name:
                    continue
                origin = _clean_name(s.get("name"))
                props = {"originName": origin} if origin else {}
                a = s.get("amount")
                if a is not None and not (isinstance(a, str) and is_unknown(a)):
                    props["amount"] = a
                if sph is not None:
                    props["solventPH"] = sph
                rels.append({"c": name, "r": rid, "t": "solvent", "p": props})
    return rels


# Adapted from CatKG_project/code/clean/relation_creat.py
def extract_article_reaction(data: Dict[str, Any]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for key, art in data.items():
        aid = key.split("_", 1)[1]
        for det in art.get("reactionDetails", []):
            rows.append({"aid": aid, "rid": str(det.get("ID"))})
    return rows


def mark_main_products(edges):
    """Mark the main products used to connect reaction steps."""
    grouped = {}
    for edge in edges:
        if edge["t"] == "product":
            grouped.setdefault(edge["r"], []).append(edge)

    def pct(value):
        if not str(value).endswith("%"):
            return 0
        try:
            return float(str(value).replace("%", ""))
        except ValueError:
            return 0

    for products in grouped.values():
        scores = [max(pct(e["p"].get("yield")), pct(e["p"].get("selectivity"))) for e in products]
        largest = max(scores)
        for edge, score in zip(products, scores):
            edge["p"]["mainProduct"] = "true" if len(products) == 1 or (largest > 0 and score == largest) else "false"
    return edges

def build_graph(cleaned, substance_map, publication_graph=None):
    data = copy.deepcopy(cleaned)
    ids = [str(rr["ID"]) for art in data.values() for rr in art["reactionDetails"]]
    _, nodes = extract_reactions(data, 1)
    # Use the same reaction IDs in the graph and literature references.
    for node, original_id in zip(nodes, ids):
        node["ID"] = original_id
    for reaction, original_id in zip((rr for art in data.values() for rr in art["reactionDetails"]), ids):
        reaction["ID"] = original_id
    articles = [article_node(int(key.split("_", 1)[1]), art["articleInformation"]) for key, art in data.items()]
    edges = mark_main_products(extract_compound(data))
    catalysts = []
    for art in data.values():
        for rr in art["reactionDetails"]:
            cat = rr["catalyst"]
            catalysts.append({"catalystName": cat.get("catalystName", "unknown"), "reactionID": str(rr["ID"]),
                **{k: v for k, v in cat.items() if k != "catalystName" and v not in (None, "unknown", [], "")}})
    counts = Counter(edge["c"] for edge in edges)
    lookup = {v.get("cid_name") or k: v for k, v in substance_map.items()}
    compounds = [{"name": name, "count": count, **lookup.get(name, {})} for name, count in sorted(counts.items())]
    journals = {}
    for article in articles:
        name = article.get("journalName")
        if name:
            journals.setdefault(name, {"journalName": name})
            if article.get("impactEffect"):
                journals[name]["impactEffect"] = article["impactEffect"]
    publication_graph = publication_graph or {}
    for journal in publication_graph.get('journal_nodes', []):
        name = journal.get('journalName')
        if name in journals:
            journals[name].update(journal)
    links = list({(x['journalName'], x['publisher']):(x) for x in publication_graph.get('journalToPublisher_relation', []) if x['journalName'] in journals}.values())
    publisher_names = {x['publisher'] for x in links}
    publishers = [x for x in publication_graph.get('publisher_nodes', []) if x['publisher'] in publisher_names]
    return {"article_nodes": articles, "reaction_nodes": nodes, "catalyst_nodes": catalysts,
            "compound_nodes": compounds, "compoundToReaction_relation": edges,
            "reactionToArticle_relation": extract_article_reaction(data),
            "journal_nodes": list(journals.values()), "publisher_nodes": publishers,
            "journalToPublisher_relation": links,
            "articleToJournal_relation": [{"aid": a["article_ID"], "journalName": a["journalName"]}
                for a in articles if a.get("journalName")]}
