#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'Normalize and flatten refinedReaction entries. Skip invalid reaction blocks and retain other fields.'

from __future__ import annotations
import json, os, re, sys, logging
from collections import OrderedDict
from typing import Any, Iterable

# Configure logging.
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
log = logging.getLogger("clean_refined")


# Schema constants.

BASE_REACTANT = {"name", "conversion"}
BASE_SOLVENT  = {"name", "amount"}
BASE_PRODUCT  = {"name", "yield", "selectivity"}
GAS_ALLOWED   = {"name", "pressure"}

THERM_COND = {"reactionType", "isHeterogeneous", "gases", "temperature", "reactionTime"}
PHOTO_COND = THERM_COND | {"lightSource", "lightWavelength", "quantumEfficiency"}
ELEC_COND  = THERM_COND | {"electrodeMaterial", "appliedPotential", "currentDensity",
                           "referenceElectrode", "counterElectrode", "energyConversionEfficiency"}

CAT_COMMON = {
    "catalystName", "catalystMixtureRatio", "promoters",
    "catalystBETSurfaceArea", "turnoverFrequency", "catalystLoading",
    "catalystSupport", "catalystShape", "catalystPoreSize",
    "catalystPoreVolume", "catPreparationMethod"
}
CAT_THERM_EXTRA = {"catAcidityAndBasicity", "catMetalOxidationStates"}

PROMOTER     = {"promoterName", "role"}
METAL_STATE  = {"metal", "oxidationState", "proportion"}

THERM_PROD = BASE_PRODUCT
PHOTO_PROD = BASE_PRODUCT
ELEC_PROD  = BASE_PRODUCT | {"faradaicEfficiency"}

SCHEMA = {
    "thermal catalysis": {"cond": THERM_COND, "prod": THERM_PROD},
    "photocatalysis"   : {"cond": PHOTO_COND,  "prod": PHOTO_PROD},
    "electrocatalysis" : {"cond": ELEC_COND ,  "prod": ELEC_PROD},
}


# Field-name mapping.


def _alts(*v: str) -> set[str]:
    return {s.lower() for s in v}

VAR: dict[str, set[str]] = {
    # Normalize aliases to name.
    "name": _alts(
        "reactant_name", "reactant name", "reactantname",
        "product_name",  "product name",  "productname",
        "gas_name",      "gas name",      "gasname",
        "solvent_name",  "solvent name",  "solventname",
        "name"
    ),

    "catalystName": _alts("catalyst_name", "catalyst name"),
    "catalystMixtureRatio": _alts("catalyst_mixture_ratio", "catalyst mixture ratio"),
    "catalystBETSurfaceArea": _alts("catalyst_betsurfacearea", "bet_surface_area", "bet surface area"),
    "turnoverFrequency": _alts("turnover_frequency", "turnover frequency"),
    "catalystLoading": _alts("catalyst_loading", "loading"),
    "catalystSupport": _alts("catalyst_support", "support"),
    "catalystShape": _alts("catalyst_shape", "shape"),
    "catalystPoreSize": _alts("catalyst_poresize", "pore_size", "pore size"),
    "catalystPoreVolume": _alts("catalyst_porevolume", "pore_volume", "pore volume"),
    "catPreparationMethod": _alts("cat_preparationmethod", "preparation_method", "preparation method"),
    "catAcidityAndBasicity": _alts("cat_acidityandbasicity", "acidity_and_basicity", "acidity and basicity"),
    "catMetalOxidationStates": _alts("cat_metaloxidationstates", "metal_oxidation_states", "metal oxidation states"),

    "reactionType": _alts("reaction_type", "reactiontype", "reaction type"),
    "isHeterogeneous": _alts("is_heterogeneous", "is heterogeneous", "isheterogeneous"),
    "temperature": _alts("temp", "reaction temperature", "temperature/°c"),
    "reactionTime": _alts("reaction_time", "reaction time", "reactiontime", "time"),
    "gases": _alts("gas", "balance_gases", "balance gases"),
    "lightSource": _alts("light_source", "light source", "lightsource"),
    "lightWavelength": _alts("light_wavelength", "light wavelength", "lightwavelength"),
    "quantumEfficiency": _alts("quantum_efficiency", "quantum efficiency", "quantumefficiency"),
    "electrodeMaterial": _alts("electrode_material", "electrode material", "electrodematerial"),
    "appliedPotential": _alts("applied_potential", "applied potential", "appliedpotential"),
    "currentDensity": _alts("current_density", "current density", "currentdensity"),
    "referenceElectrode": _alts("reference_electrode", "reference electrode", "referenceelectrode"),
    "counterElectrode": _alts("counter_electrode", "counter electrode", "counterelectrode"),
    "energyConversionEfficiency": _alts("energy_conversion_efficiency", "energy conversion efficiency",
                                        "energyconversionefficiency"),

    "selectivity": _alts("sel", "select", "selectivity (%)", "sel%"),
    "yield": _alts("yield (%)", "yield%"),
    "faradaicEfficiency": _alts("faradaic_efficiency", "faradaic efficiency (%)", "fe"),

    "promoterName": _alts("promoter_name", "promoter", "name"),
    "role": _alts("role"),
    "metal": _alts("metal"),
    "oxidationState": _alts("oxidation_state", "oxidationstate"),
    "proportion": _alts("proportion"),
}

MAP: dict[str, str] = {alt: std for std, alts in VAR.items() for alt in alts}
MAP.update({k.lower(): k for k in VAR})

def std(k: str) -> str:  # Normalize a field name.
    return MAP.get(k.strip().lower(), k.strip())


def keep(d: Any, allowed: set[str]) -> dict | None:
    'Keep allowed fields; return None for an invalid dictionary block.'
    if not isinstance(d, dict):
        return None
    return {k: d.get(k, "unknown") for k in allowed}


def _append_clean(out: list[dict], item: dict | None, what: str):
    'Store a cleaned item; log and skip invalid entries.'
    if item is None:
        log.warning("Skip dirty %s entry", what)
    else:
        out.append(item)


def scrub_list(lst: Iterable[Any] | None, allowed: set[str]) -> list[dict]:
    'Normalize list entries and skip unsupported types.'
    out: list[dict] = []
    for item in lst or []:
        # Normalize a plain string entry.
        if isinstance(item, str):
            if allowed is GAS_ALLOWED:
                _append_clean(out, {"name": item, "pressure": "unknown"}, "gas")
                continue
            if allowed is BASE_SOLVENT:
                _append_clean(out, {"name": item, "amount": "unknown"}, "solvent")
                continue
            if allowed is PROMOTER:
                _append_clean(out, {"promoterName": item, "role": "unknown"}, "promoter")
                continue
            if allowed is METAL_STATE:
                _append_clean(out, {"metal": item, "oxidationState": "unknown", "proportion": "unknown"}, "metalState")
                continue
            log.warning("Skip dirty list element (str) → %r", item)
            continue

        if isinstance(item, dict):
            cleaned = {std(k): v for k, v in item.items()}
            _append_clean(out, keep(cleaned, allowed), "dict")
        else:
            log.warning("Skip dirty list element (type=%s)", type(item).__name__)
    return out


def _clean_one(rr: dict) -> dict | None:
    'Normalize one refinedReaction entry; return None for an invalid block.'
    try:
        rr = {std(k): v for k, v in rr.items()}
        cond = rr.get("conditions", {}) or {}
        rtype = str(cond.get("reactionType", "")).lower()
        if rtype not in SCHEMA:
            if any(k in cond for k in ("lightSource", "quantumEfficiency")):
                rtype = "photocatalysis"
            elif any(k in cond for k in ("electrodeMaterial", "appliedPotential", "currentDensity")):
                rtype = "electrocatalysis"
            else:
                rtype = "thermal catalysis"
        spec = SCHEMA[rtype]

        out: OrderedDict[str, Any] = OrderedDict()
        # reactants
        reac = rr.get("reactants", {}) or {}
        out["reactants"] = {
            "reactantList": scrub_list(reac.get("reactantList"), BASE_REACTANT),
            "reactantMixtureRatio": reac.get("reactantMixtureRatio", "unknown"),
        }
        # catalyst
        allowed_cat = CAT_COMMON | (CAT_THERM_EXTRA if rtype == "thermal catalysis" else set())
        cat = rr.get("catalyst", {})
        cat_blk = keep(cat, allowed_cat)
        if cat_blk is None:  # Skip the reaction when the catalyst block has an invalid type.
            raise ValueError("catalyst field malformed")
        cat_blk["promoters"] = scrub_list(cat.get("promoters"), PROMOTER)
        if "catMetalOxidationStates" in allowed_cat:
            cat_blk["catMetalOxidationStates"] = scrub_list(cat.get("catMetalOxidationStates"), METAL_STATE)
        out["catalyst"] = cat_blk
        # conditions
        gases = scrub_list(cond.get("gases"), GAS_ALLOWED)
        cond_blk = keep(cond, spec["cond"] - {"gases"})
        if cond_blk is None:
            raise ValueError("conditions block malformed")
        cond_blk["gases"] = gases
        out["conditions"] = cond_blk
        # solvents
        solv = rr.get("solvents", {}) or {}
        out["solvents"] = {
            "solventList": scrub_list(solv.get("solventList"), BASE_SOLVENT),
            "solventPH": solv.get("solventPH", "unknown"),
        }
        # products
        out["products"] = scrub_list(rr.get("products"), spec["prod"])
        return out
    except Exception as e:
        log.warning("Skip refinedReaction due to error: %s", e)
        return None


def clean_refined(rr: Any) -> list[dict] | dict | None:
    'Accept supported list or dictionary inputs.'
    if isinstance(rr, list):
        cleaned_list = [_clean_one(x) for x in rr if isinstance(x, dict)]
        return [x for x in cleaned_list if x is not None]
    if isinstance(rr, dict):
        return _clean_one(rr)
    return None


REG = re.compile(r"refined[\s_]?reaction", re.I)

KEY_ALIASES = {
    "article_information": "articleInformation",
    "article info": "articleInformation",
    "articleinfo": "articleInformation",
    "reaction_details": "reactionDetails",
    "reactiondetails": "reactionDetails",
}


def _normalize_keys(node: dict) -> OrderedDict:
    out: OrderedDict[str, Any] = OrderedDict()
    for k, v in node.items():
        if isinstance(k, str):
            kl = k.strip().lower()
            mapped = KEY_ALIASES.get(kl)
            if mapped and mapped not in node:
                out[mapped] = v
                continue
        out[k] = v
    return out


def _normalize_reaction_details(node: OrderedDict) -> OrderedDict:
    rd = node.get("reactionDetails")
    if isinstance(rd, dict) and all(str(k).isdigit() for k in rd.keys()):
        items = []
        for idstr, blk in rd.items():
            if isinstance(blk, dict):
                items.append({"ID": int(idstr), **blk})
            else:
                items.append({"ID": int(idstr), "value": blk})
        node["reactionDetails"] = items
    return node

def walk(node: Any):
    'Recursively normalize, flatten, and filter refinedReaction records.'
    if isinstance(node, dict):
        node = _normalize_reaction_details(_normalize_keys(node))
        # Flatten cleaned refinedReaction entries at this level.
        for k, v in list(node.items()):
            if REG.fullmatch(k):
                items = clean_refined(v)
                if not items:
                    # Drop the parent when no refined reaction survives.
                    log.warning("Node dropped: all refinedReaction invalid")
                    return None
                items_list = items if isinstance(items, list) else [items]
                result = []
                for cleaned in items_list:
                    new_node: OrderedDict[str, Any] = OrderedDict()
                    # Retain the other fields.
                    for k2, v2 in node.items():
                        if k2 == k:
                            continue
                        w = walk(v2)
                        if w is not None:
                            new_node[k2] = w
                    # Merge the cleaned fields.
                    new_node.update(cleaned)
                    result.append(new_node)
                return result
        # Recurse into dictionary values.
        res: OrderedDict[str, Any] = OrderedDict()
        for k, v in node.items():
            w = walk(v)
            if w is not None:
                res[k] = w
        return res

    if isinstance(node, list):
        out = []
        for item in node:
            w = walk(item)
            if isinstance(w, list):
                out.extend(w)
            elif w is not None:
                out.append(w)
        return out

    # Return scalar values unchanged.
    return node


def main() -> None:
    files = [a for a in sys.argv[1:] if not a.startswith("-")]
    if not files:
        print("Usage: normalize_structure.py <input.json> [input.json ...]")
        sys.exit(1)
    for fp in files:
        with open(fp, encoding="utf-8") as f:
            raw = json.load(f)
        cleaned = walk(raw)
        if cleaned is None:
            log.error("%s → nothing left after cleaning; skip output", fp)
            continue
        out_path = os.path.splitext(fp)[0] + "_structured.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(cleaned, f, ensure_ascii=False, indent=2)
        print("[OK] wrote", out_path)

if __name__ == "__main__":
    main()
