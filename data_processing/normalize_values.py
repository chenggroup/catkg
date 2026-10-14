#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'Normalize performance, BET, pressure, temperature, and time fields in flattened reactionDetails records. Skip invalid blocks and log changes.'

from __future__ import annotations
import argparse
import json, os, re, sys, logging, traceback
from pathlib import Path
from typing import Any


LOG = logging.getLogger("clean-val")
LOG.setLevel(logging.INFO)
LOG.addHandler(logging.StreamHandler(sys.stdout))   # Send log records to stdout.


SUBSTANCE_MAP: dict[str, dict[str, Any]] = {}


def load_substance_map(path: Path) -> None:
    global SUBSTANCE_MAP
    data = json.load(path.open(encoding="utf-8"))
    SUBSTANCE_MAP = {
        k: {"cid": v.get("cid"), "cidName": v.get("cid_name"), "chebiID": v.get("chebi_id")}
        for k, v in data.items() if isinstance(v, dict)
    }


UNAVAIL = {"unknown", "n/a", "not provided", "not mentioned", "not specified", "unspecified"}

def _name_ok(n: str) -> bool:
    return bool(n and n.strip().lower() not in UNAVAIL)


def _log_set(obj: dict, key: str, new: str, where: str):
    old = obj.get(key, "⟂None")
    if old != new:
        LOG.info("%s.%s : '%s' -> '%s'", where, key, old, new)
        obj[key] = new


def _add_ids(d: dict, where: str):
    nm = d.get("name", "")
    if _name_ok(nm) and (m := SUBSTANCE_MAP.get(nm.strip())):
        for k, v in m.items():
            if v and k not in d:
                LOG.info("%s.%s : <insert> '%s'", where, k, v)
                d[k] = v


_pct_pat = [
    r'(\d+(?:\.\d+)?)[–-](\d+(?:\.\d+)?)%',
    r'(\d+(?:\.\d+)?)[ ]*±[ ]*\d+(?:\.\d+)?%',
    r'(\d+(?:\.\d+)?)%'
]
_last = lambda s: next((m.group(len(m.groups())//2) for p in _pct_pat if (m := re.search(p, s))), None)


def _pct(x: Any) -> str:
    s = str(x).strip().replace('\u200b', '')
    s = re.sub('%+', '%', s)
    if "%" in s:
        s = re.sub(r'\s*%\s*', '%', s).replace('>', '')
        ret = "unknown" if '<' in s else (_last(s) + "%" if _last(s) else "unknown")
    else:
        try:
            f = float(s)
            if 0 <= f <= 1:
                ret = f"{f*100:.2f}%"
            elif 1 < f <= 100:
                ret = f"{f}%"
            else:
                ret = "unknown"
        except Exception:
            sl = s.lower()
            if sl.startswith(("full", "quantitative")):
                ret = "100%"
            elif sl.startswith("trace"):
                ret = "0%"
            else:
                ret = "unknown"
    return re.sub('%+', '%', ret)

conv = yield_ = sel = _pct

# Normalize BET units.
BET_SUB = re.compile('|'.join(map(re.escape, [
    "m²/g", "m2/g", "m2g-1", "m2·g−1", "m²g⁻¹", "m2·g–1", "m2·g-1", "m²g⁻1"
])))

def _bet(v: Any) -> str:
    if not v or str(v).lower() == "unknown":
        return "unknown"
    s = BET_SUB.sub('m2·g−1', str(v))
    for p in [
        r'(\d+(?:\.\d+)?)[–-](\d+(?:\.\d+)?)',
        r'(\d+(?:\.\d+)?)[ ]*±[ ]*\d',
        r'(\d+(?:\.\d+)?)'
    ]:
        if (m := re.search(p, s)):
            return f"{m.group(len(m.groups())//2)} m2·g−1"
    return "unknown"

# Normalize pressure units.
UNIT = {
    'bar': 1, 'bars': 1, 'barg': 1, 'barG': 1,
    'mbar': 0.001, 'mBar': 0.001,
    'atm': 1.01325, 'atmosphere': 1.01325,
    'kPa': 0.01, 'MPa': 10, 'GPa': 1e4, 'Pa': 1e-5,
    'Torr': 0.00133322, 'torr': 0.00133322,
    'mTorr': 1.33322e-6,
    'psi': 0.0689476, 'psig': 0.0689476,
    'MPag': 10, 'mmHg': 1.01325/760, 'mm Hg': 1.01325/760
}
ATM_STR = {s.lower() for s in [
    "ambient pressure", "ambient", "atmospheric pressure", "atmospheric",
    "ordinary pressure", "normal pressure", "ambient atmosphere",
    "atmosphere pressure", "1 atmosphere", "ambient air pressure",
    "room atmosphere", "nearly atmospheric pressure",
    "balance to atmospheric pressure", "ambient (balloon)",
    "close to atmospheric pressure", "atmospheric flow",
    "near ambient", "standard pressure", "nearly atmospheric", "room pressure"
]}

def _fmt_bar(v: float) -> str:
    return f"{v:.2f} bar" if v >= 1 else f"{v:.3g} bar"


def _clean_pressure(val: Any) -> str:
    if not val:
        return "unknown"
    s = str(val).strip()
    if s.lower() in ATM_STR or s == "760 mm Hg":
        return "1.01325 bar"
    sci = re.search(r'([\d\.]+(?:[Ee][+\-]?\d+| × 10\^?\-?\d+))\s*([a-zA-Z]+)', s)
    if sci and (unit := sci.group(2)) in UNIT:
        num = eval(sci.group(1).replace('×', '*').replace('^', '**'))
        return _fmt_bar(num * UNIT[unit])
    rng = re.search(r'(\d+(?:\.\d+)?)\s*(?:±|to|-|–|≤|≥)\s*(\d+(?:\.\d+)?)\s*([a-zA-Z]+)', s)
    if rng and (unit := rng.group(3)) in UNIT:
        num = max(float(rng.group(1)), float(rng.group(2)))
        return _fmt_bar(num * UNIT[unit])
    mv = re.search(r'(\d+(?:\.\d+)?)\s*([a-zA-Z]+)', s)
    if mv and (unit := mv.group(2)) in UNIT:
        return _fmt_bar(float(mv.group(1)) * UNIT[unit])
    return "unknown"

# Other value cleaners.

def _bval(v: Any) -> str:
    t = str(v).strip().lower()
    if t in {"true", "yes", "heterogeneous"}:
        return "true"
    if t in {"false", "no", "homogeneous"}:
        return "false"
    return "unknown"


def _temp(v: Any) -> str:
    s = str(v).strip().lower()
    if s in {"room temperature", "rt", "ambient", "ambient temperature", "298 k", "r.t."}:
        return "298 K"
    if s in {"0 °c", "0°c", "0"}:
        return "273 K"
    rules = [
        (r'(\d+(?:\.\d+)?)[–-](\d+(?:\.\d+)?)\s*k', lambda m: f"{m.group(2)} K"),
        (r'(\d+(?:\.\d+)?)[–-](\d+(?:\.\d+)?)\s*°c', lambda m: f"{float(m.group(2))+273:.0f} K"),
        (r'([-+]?\d+(?:\.\d+)?)\s*°c', lambda m: f"{float(m.group(1))+273:.0f} K"),
        (r'([-+]?\d+(?:\.\d+)?)\s*k', lambda m: m.group(1)+" K"),
    ]
    for pat, fn in rules:
        if (m := re.search(pat, s)):
            return fn(m)
    try:
        return f"{float(s)+273:.0f} K"
    except Exception:
        return "unknown"


def _rtime(v: Any) -> str:
    s = str(v).strip().lower()
    if s in {"unknown", "until complete", "appropriate time"}:
        return "unknown"
    if s == "overnight":
        return "12 h"
    m = re.search(r'(\d+(?:\.\d+)?)[ ]*(ms|μs|us|ps|s|sec|secs|seconds|min|mins|minutes|h|hr|hrs|hours|d|day|days)', s)
    if not m:
        return "unknown"
    num, unit = float(m.group(1)), m.group(2)
    if unit.startswith('min'):
        return f"{num/60:.2f} h"
    if unit.startswith('h'):
        return f"{num} h"
    if unit.startswith('d'):
        return f"{num*24} h"
    if unit.startswith('s'):
        return f"{num} s"
    return "unknown"

VALID_RT = {"thermal catalysis", "photocatalysis", "electrocatalysis"}

# Clean one reaction.

def _clean_rr(rr: dict, where: str):
    # reactants
    rlist = rr.get("reactants", {}).get("reactantList", [])
    for i, r in enumerate(rlist, 1):
        if not isinstance(r, dict):
            LOG.warning("Skip non-dict reactant at %s.reactants[%d]", where, i)
            continue
        _log_set(r, "conversion", conv(r.get("conversion", "unknown")), f"{where}.reactants[{i}]")
        _add_ids(r, f"{where}.reactants[{i}]")

    # products
    for i, pd in enumerate(rr.get("products", []), 1):
        if not isinstance(pd, dict):
            LOG.warning("Skip non-dict product at %s.products[%d]", where, i)
            continue
        _log_set(pd, "yield",       yield_(pd.get("yield", "unknown")),       f"{where}.products[{i}]")
        _log_set(pd, "selectivity", sel(pd.get("selectivity", "unknown")),   f"{where}.products[{i}]")
        _add_ids(pd, f"{where}.products[{i}]")

    # gases
    for i, g in enumerate(rr.get("conditions", {}).get("gases", []), 1):
        if not isinstance(g, dict):
            LOG.warning("Skip non-dict gas at %s.conditions.gases[%d]", where, i)
            continue
        _log_set(g, "pressure", _clean_pressure(g.get("pressure", "unknown")), f"{where}.conditions.gases[{i}]")
        _add_ids(g, f"{where}.conditions.gases[{i}]")

    # solvents
    for i, s in enumerate(rr.get("solvents", {}).get("solventList", []), 1):
        if not isinstance(s, dict):
            LOG.warning("Skip non-dict solvent at %s.solvents[%d]", where, i)
            continue
        _add_ids(s, f"{where}.solvents[{i}]")

    # BET & conditions
    cat = rr.get("catalyst", {}) or {}
    if isinstance(cat, dict):
        _log_set(cat, "catalystBETSurfaceArea", _bet(cat.get("catalystBETSurfaceArea", "unknown")), f"{where}.catalyst")

    cond = rr.get("conditions", {}) or {}
    if isinstance(cond, dict):
        pc = f"{where}.conditions"
        _log_set(cond, "isHeterogeneous", _bval(cond.get("isHeterogeneous", "unknown")), pc)
        _log_set(cond, "temperature", _temp(cond.get("temperature", "unknown")), pc)
        _log_set(cond, "reactionTime", _rtime(cond.get("reactionTime", "unknown")), pc)
        rt = str(cond.get("reactionType", "unknown")).strip().lower()
        _log_set(cond, "reactionType", rt if rt in VALID_RT else "others", pc)

# Traverse the JSON input.

def _traverse(node: Any, p: str = "$"):
    if isinstance(node, dict):
        if {"reactants", "products", "conditions"}.issubset(node):
            _clean_rr(node, p)
        for k, v in node.items():
            _traverse(v, f"{p}.{k}")
    elif isinstance(node, list):
        for idx, it in enumerate(node, 1):
            _traverse(it, f"{p}[{idx}]")


def main() -> None:
    try:
        ap = argparse.ArgumentParser(description="Clean numeric values in structured JSON")
        ap.add_argument("input_json", help="Input JSON file")
        ap.add_argument(
            "--substance-map",
            default=str(Path(__file__).resolve().parents[3] / "data/compound_names.json"),
            help="Substance map JSON path",
        )
        args = ap.parse_args()

        infile = Path(args.input_json)
        map_path = Path(args.substance_map).expanduser().resolve()
        if not map_path.exists():
            raise FileNotFoundError(map_path)
        load_substance_map(map_path)

        data = json.loads(infile.read_text(encoding="utf-8"))
        _traverse(data)
        out = infile.with_name(infile.stem + "_cleaned.json")
        out.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
        print("[OK] wrote", out)
    except Exception:
        LOG.error("normalize_values.py error:\n%s", traceback.format_exc())
        sys.exit(2)

if __name__ == "__main__":
    main()
