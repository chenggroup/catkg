"""Score individual reactions and paths, including their A/B/C contributions."""
import ast
import re
from datetime import datetime

def calculate_score_impact_effect(impact_effect: float) -> float:
    max_impact_effect = 50.5
    return round((impact_effect / max_impact_effect) * 200, 2)

def calculate_score_date(date_str: str) -> float:
    start_date = datetime(2000, 1, 1)
    end_date   = datetime(2024, 7, 30)
    for fmt in ("%Y/%m/%d","%Y-%m-%d","%Y/%m","%Y-%m"):
        try:
            date = datetime.strptime(date_str, fmt)
            break
        except ValueError:
            continue
    else:
        return 0.0
    if date >= end_date:
        return 100.0
    if date < start_date:
        return 0.0
    total_days   = (end_date - start_date).days
    days_from_end = (end_date - date).days
    return round((1 - days_from_end / total_days) * 100, 2)

def score_citations(citations) -> float:
    try:
        c = int(citations)
    except:
        c = 0
    return 200.0 if c > 100 else max(c, 0)

def calculate_catalyst_score(c: dict) -> int:
    keys = [
        "catalystLoading", "catalystSupport", "catPreparationMethod",
        "promoters", "catMetalOxidationStates", "catalystBETSurfaceArea",
        "turnoverFrequency", "catalystShape", "catalystPoreSize", "catalystPoreVolume"
    ]
    score_each = 20
    return sum(score_each for k in keys if k in c and c[k] not in (None, "unknown", [], ""))

def calculate_conditions_score(cond: dict) -> int:
    return 100 if cond.get("temperature") not in (None, "unknown", "") else 0

def check_gases_conditions(cond: dict) -> int:
    gs = cond.get("gases")
    if isinstance(gs, list) and gs:
        valid = [g for g in gs if g.get("name") not in (None, "unknown", "") and g.get("pressure") not in (None, "unknown", "")]
        if len(valid) == len(gs):
            return 100
    return 0

def calculate_reaction_time_score(cond: dict) -> int:
    return 50 if cond.get("reactionTime") not in (None, "unknown", "") else 0

def _pct(value) -> float:
    if value in (None, "", "unknown"):
        return 0.0
    if isinstance(value, (int, float)):
        return float(value)/100 if value > 1 else float(value)
    m = re.search(r"[-+]?\d*\.?\d+", str(value))
    if m:
        return float(m.group())/100
    return 0.0

def calculate_reaction_performance_score(rxn: dict) -> float:
    reactants = rxn.get("reactants", {}).get("reactantList", [])
    products  = rxn.get("products", [])
    max_conv = max((_pct(r.get("conversion")) for r in reactants), default=0.0)
    total = 0.0
    for p in products:
        yld = _pct(p.get("yield"))
        sel = _pct(p.get("selectivity"))
        total += 400 * max(sel * max_conv, yld) + 100 * (2 * sel + max_conv + yld)
    return total

def add_scores_and_single_flag(data: dict) -> None:
    for article in data.values():
        details = article.get("reactionDetails", [])
        single_flag = len(details) == 1
        info = article.get("articleInformation", {})
        impact_score   = calculate_score_impact_effect(float(info.get("impactEffect", 0)))
        date_score     = calculate_score_date(str(info.get("pubDate", "2000-01-01")))
        citation_score = score_citations(info.get("citations", 0))
        for rxn in details:
            cond = rxn.get("conditions", {})
            cond["singleReaction"] = "true" if single_flag else "false"
            rxn["conditions"] = cond
            c_score = calculate_catalyst_score(rxn.get("catalyst", {}))
            cond_score = calculate_conditions_score(cond)
            time_score = calculate_reaction_time_score(cond)
            gases_score = check_gases_conditions(cond)
            perf_score = calculate_reaction_performance_score(rxn)
            total = impact_score + date_score + citation_score + c_score + cond_score + time_score + gases_score + perf_score
            for p in rxn.get("products", []):
                p["score"] = f"{total:.2f}"

CATALYST_FIELDS = ["catalystLoading", "catalystSupport", "catPreparationMethod", "promoters",
                  "catMetalOxidationStates", "catalystBETSurfaceArea", "turnoverFrequency",
                  "catalystShape", "catalystPoreSize", "catalystPoreVolume"]
SIMPLE = {"CO", "H2", "syngas", "H2O", "CO2", "Carbon Dioxide", "Carbon Monoxide",
          "Carbon monoxide", "Carbon dioxide", "Hydrogen", "Oxygen", "Water",
          "hydrogen", "O2", "oxygen", "Nitrogen", "water", "N2"}


def as_list(value):
    return ast.literal_eval(value) if isinstance(value, str) else (value or [])


def single_breakdown(info, reaction):
    cond, cat = reaction.get("conditions", {}), reaction.get("catalyst", {})
    terms = {
        "A1_journal_impact": calculate_score_impact_effect(float(info.get("impactEffect", 0) or 0)),
        "A2_publication_date": calculate_score_date(str(info.get("pubDate", "2000-01-01"))),
        "A3_citations": score_citations(info.get("citations", 0)),
        "A4_catalyst_completeness": calculate_catalyst_score(cat),
        "A5_temperature_available": calculate_conditions_score(cond),
        "A6_all_gas_pressures_available": check_gases_conditions(cond),
        "A7_reaction_time_available": calculate_reaction_time_score(cond),
        "A8_performance": calculate_reaction_performance_score(reaction),
    }
    max_conversion = max((_pct(r.get("conversion")) for r in reaction.get("reactants", {}).get("reactantList", [])), default=0)
    products = []
    for product in reaction.get("products", []):
        y, s = _pct(product.get("yield")), _pct(product.get("selectivity"))
        products.append({"name": product["name"], "yield_cleaned": product.get("yield"),
            "selectivity_cleaned": product.get("selectivity"), "yield_fraction": y,
            "selectivity_fraction": s, "max_conversion_fraction": max_conversion,
            "contribution": 400 * max(s * max_conversion, y) + 100 * (2*s + max_conversion + y)})
    return {"terms": terms, "products": products,
            "known_catalyst_fields": [k for k in CATALYST_FIELDS if cat.get(k) not in (None, "unknown", [], "")],
            "total_unrounded": sum(terms.values()), "total": round(sum(terms.values()), 2)}


def pressure(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    s = str(value).strip().lower()
    if s in ("", "unknown"):
        return None
    for unit in [" bar", "bar", " kpa", "mpa"]:
        if s.endswith(unit):
            s = s[:-len(unit)].strip()
            break
    try:
        return float(s)
    except ValueError:
        return None


def pair_breakdown(a, b):
    out = {f"B{i}": 0.0 for i in range(1, 9)}
    for number, field, weight in [(1, "isHeterogeneous", 50), (2, "catalystName", 200), (3, "catalystSupport", 100)]:
        av, bv = a.get(field, "").lower(), b.get(field, "").lower()
        if av not in ("unknown", "") and bv not in ("unknown", "") and av == bv:
            out[f"B{number}"] = weight
    promoters = lambda r: {p["promoterName"].lower() for p in as_list(r.get("promoters", "[]")) if "promoterName" in p}
    pa, pb = promoters(a), promoters(b)
    if pa and pb and pa == pb:
        out["B4"] = 100
    ta, tb = a.get("temperature", "unknown").replace(" K", ""), b.get("temperature", "unknown").replace(" K", "")
    temperature_difference = None
    try:
        if ta not in ("unknown", "") and tb not in ("unknown", ""):
            temperature_difference = abs(float(ta)-float(tb))
            out["B5"] = 30000/(temperature_difference + 30)
    except ValueError:
        pass
    ga = [g for g in as_list(a.get("gases", [])) if isinstance(g, dict) and g.get("name", "").strip().lower() not in ("", "unknown")]
    gb = [g for g in as_list(b.get("gases", [])) if isinstance(g, dict) and g.get("name", "").strip().lower() not in ("", "unknown")]
    common = {g["name"].lower() for g in ga} & {g["name"].lower() for g in gb}
    if common:
        out["B6"] = 100
        for name in common:
            p1 = pressure(next(g for g in ga if g["name"].lower() == name).get("pressure"))
            p2 = pressure(next(g for g in gb if g["name"].lower() == name).get("pressure"))
            if p1 is not None and p2 is not None and p1+p2 > 0:
                out["B7"] = max(out["B7"], 200*(1-abs(p1-p2)/(p1+p2)))
    if set(as_list(a.get("solventList", []))) & set(as_list(b.get("solventList", []))):
        out["B8"] = 100
    return {"terms": out, "temperature_difference_K": temperature_difference, "total": sum(out.values())}


def path_breakdown(result):
    """Calculate A/B/C for a path listed in backward retrieval order."""
    steps = result["reaction"]
    pairs = [pair_breakdown(a, b) for a, b in zip(steps, steps[1:])]
    compounds = set()
    for step in steps:
        compounds.update(as_list(step.get("reactantList", [])))
        compounds.update(g.get("name", "unknown") for g in as_list(step.get("gases", [])) if isinstance(g, dict))
        compounds.update(as_list(step.get("solventList", [])))
    c1 = 500 if all(r in SIMPLE for r in as_list(steps[-1].get("reactantList", []))) else 0
    c2 = 10*(15-len(compounds)) if len(compounds) <= 15 else 0
    a = sum(float(s) if s is not None else 0 for s in result["score"])
    b, c = sum(p["total"] for p in pairs), c1+c2
    # Accumulate in the scoring order: adjacent-step terms, reaction scores, C1/C2.
    accumulated = 0.0
    for pair in pairs:
        for value in pair["terms"].values():
            accumulated += value
    for value in result["score"]:
        accumulated += float(value) if value is not None else 0.0
    accumulated += c1
    accumulated += c2
    total = round(accumulated, 2)
    return {"A": a, "A_step_scores_reverse_order": result["score"], "B": b, "B_pairs": pairs,
            "C": c, "C1_simple_feed": c1, "C2_fewer_compounds": c2,
            "C2_counted_compounds": sorted(compounds), "total": total,
            "excluded_by_H2_O2_rule": "H2" in compounds and "O2" in compounds}


def filter_data(data):
    filtered_data = []

    for item in data:
        # Collect every compound name in the path.
        path_items = item.get('path',[])
        substance_names = []


        for path_item in path_items:

            if isinstance(path_item, list):
                substance_names.extend(path_item)
            else:
                substance_names.append(path_item)

        # Exclude paths with repeated names.
        if len(substance_names) == len(set(substance_names)):
            filtered_data.append(item)

    return filtered_data

def calculate_temperature_difference(filtered_data, max_results=100):
    'Sort candidate paths by absolute temperature difference and retain at most max_results.'
    path_with_temp_diff = []
    for item in filtered_data:
        temps = [float(reaction.get('temperature').replace(" K", "")) for reaction in item.get('reaction', []) if reaction.get('temperature')]
        if temps:
            temp_diff = max(temps) - min(temps)
            path_with_temp_diff.append((temp_diff, item))
    
    path_with_temp_diff.sort(key=lambda x: x[0])


    # Retain the requested number of pairs with the smallest temperature differences.
    filtered_paths = [path for _, path in path_with_temp_diff[:max_results]]
    print("Temperature-shortlisted paths:",len(filtered_paths))
    return filtered_paths

def calculate_score_between_reactions(candidates):
    """Rank candidates by A+B+C using the same terms as the score explanation."""
    ranked = []
    for candidate in candidates:
        detail = path_breakdown(candidate)
        if not detail["excluded_by_H2_O2_rule"]:
            ranked.append((candidate, detail["total"]))
    return sorted(ranked, key=lambda pair: pair[1], reverse=True)
