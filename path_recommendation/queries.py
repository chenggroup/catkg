"""Build backward and branch-and-bound Cypher queries for pathway retrieval.

General queries require APOC for graph properties stored as JSON lists.
The fixed 40-reaction sample queries below use parameters and require no APOC.
"""
import os
import json
from pathlib import Path
from typing import List, Dict, Any, Optional

substance_map = {}

def get_substance_name(name, substance_map_path=None):
    """Return the canonical compound name, its recorded alias, or None if absent."""
    global substance_map
    if substance_map_path is None:
        substance_map_path = Path(__file__).resolve().parents[1] / "data/compound_names.json"
    # Load the bundled compound map.
    if not os.path.exists(substance_map_path):
        raise FileNotFoundError(f"Substance map file not found: {substance_map_path}")

    if not substance_map:
        with open(substance_map_path, 'r', encoding='utf-8') as f:
            json_data = json.load(f)
            for item in json_data.values():
                cid_name = item.get('cid_name', '')
                if cid_name and cid_name.strip() and cid_name not in substance_map:
                    substance_map[cid_name] = item
            for key, value in json_data.items():
                if key not in substance_map:
                    substance_map[key] = value
    
    if not name or not name.strip():
        return None

    substance_info = substance_map.get(name)
    if substance_info:
        cid_name = substance_info.get('cid_name')
        if cid_name and cid_name.strip():
            return cid_name
        else:
            return name
    else:
        return None

def generate_optimized_cypher_query(target_compound, step_conditions, path_limitation, mainProduct="true", to_jlist = True):
    def jlist(field):
        if not to_jlist:
            return field
        return f"apoc.convert.fromJsonList({field})"

    query = f"""
    MATCH (s:Compound {{name: '{get_substance_name(target_compound)}'}})
    MATCH (s)<-[r:has_product]-(m:Reaction)
    MATCH (m)<-[:is_article_of]-(:Article)-[:is_published_in]->(j:Journal)
    WHERE m.temperature IS NOT NULL 
    AND NOT '{get_substance_name(target_compound)}' IN {jlist("m.reactantList")}
    """
    if mainProduct == "true":
        query += "AND r.mainProduct = 'true' "


    conditions = step_conditions[0]

    if "yield_limitation" in conditions and conditions["yield_limitation"] > 0:
        query += f"AND r.yield IS NOT NULL AND toFloat(REPLACE(r.yield, '%', '')) > {conditions['yield_limitation']} "
    if "selectivity_limitation" in conditions and conditions["selectivity_limitation"] > 0:
        query += f"AND r.selectivity IS NOT NULL AND toFloat(REPLACE(r.selectivity, '%', '')) > {conditions['selectivity_limitation']} "
    if "score" in conditions and conditions["score"] > 0:
        query += f"AND r.score IS NOT NULL AND toFloat(r.score) > {conditions['score']} "
    if conditions.get("isHeterogeneous", "") == "true":
        query += f"AND m.isHeterogeneous = 'true' "

    if "include_ids" in conditions and conditions["include_ids"]:
        for id_val in conditions['include_ids']:
            query += f"AND m.ID = '{id_val}' "
    if "exclude_ids" in conditions and conditions["exclude_ids"]:
        for id_val in conditions["exclude_ids"]:
            query += f"AND m.ID <> '{id_val}' "
    if "must_contain" in conditions:
        for compound in conditions["must_contain"]:
            query += f"AND '{get_substance_name(compound)}' IN {jlist('m.reactantList')} "

    if "must_contain_products" in conditions:
        for compound in conditions["must_contain_products"]:
            query += f"AND '{get_substance_name(compound)}' IN {jlist('m.productList')} "

    if "must_not_contain" in conditions:
        for compound in conditions["must_not_contain"]:
            std_name = get_substance_name(compound)
            query += f"""
            AND NOT '{std_name}' IN {jlist("m.reactantList")}
            AND NOT '{std_name}' IN {jlist("m.solventList")}
            AND ALL(g IN apoc.convert.fromJsonList(m.gases) WHERE g.name IS NULL OR g.name <> '{std_name}')
            """

    if "exclude_products" in conditions:
        for compound in conditions["exclude_products"]:
            query += f"AND NOT '{get_substance_name(compound)}' IN {jlist('m.productList')} "

    if "max_product_names" in conditions and conditions['max_product_names'] not in [False, 0]:
        query += f"AND size({jlist('m.productList')}) < {conditions['max_product_names']+1} "
    if "max_reactant_names" in conditions and conditions['max_reactant_names'] not in [False, 0]:
        query += f"AND size({jlist('m.reactantList')}) < {conditions['max_reactant_names']+1} "
    if conditions.get("isSingleReaction", "") == "true":
        query += f"AND m.singleReaction = 'true' "
    if "temperature" in conditions and conditions["temperature"] != False:
        temp_min, temp_max = conditions["temperature"]
        query += f"""
        AND toFloat(REPLACE(m.temperature, ' K', '')) >= {temp_min}
        AND toFloat(REPLACE(m.temperature, ' K', '')) <= {temp_max}
        """
    if "citations" in conditions and conditions["citations"] > 0:
        query += f"AND toInteger(m.citations) >= {conditions['citations']} "
    if "impactEffect" in conditions and conditions["impactEffect"] > 0:
        query += f"AND toFloat(m.impactEffect) >= {conditions['impactEffect']} "
    if "reactionType" in conditions and conditions["reactionType"]:
        query += f"AND m.reactionType = '{conditions['reactionType']}' "
    if "catalyst" in conditions and conditions["catalyst"]:
        query += f"AND m.catalystName CONTAINS '{conditions['catalyst']}' "
    if "pubDate" in conditions and conditions["pubDate"]:
        query += f"AND toInteger(substring(m.pubDate, 0, 4)) >= {conditions['pubDate']} "

    for step_index, cond in enumerate(step_conditions):
        if step_index == 0:
            continue
        previous_step_label = f"m{step_index - 1}" if step_index > 1 else "m"
        step_label = f"m{step_index}"
        reactant_rel = f"r{step_index*2 - 1}"
        product_rel = f"r{step_index*2}"

        query += f"""
        MATCH ({previous_step_label})<-[{reactant_rel}:is_reactant_of]-(s{step_index}:Compound)
        MATCH (s{step_index})<-[{product_rel}:has_product]-({step_label}:Reaction)
        MATCH ({step_label})<-[:is_article_of]-(:Article)-[:is_published_in]->(j{step_index}:Journal)
        WHERE {step_label}.temperature IS NOT NULL 
        AND NOT '{get_substance_name(target_compound)}' IN {jlist(f"{step_label}.reactantList")}        
        """


        if mainProduct == "true":
            query += f"AND {product_rel}.mainProduct = 'true' "

        # Apply the bridge-compound constraint.
        # Read must_contain from the next step toward the target.
        # Require the current upstream product to belong to that list.
        # The bridge must be a main product regardless of the global flag.
        prev_conditions = step_conditions[step_index - 1]
        if "must_contain" in prev_conditions and prev_conditions["must_contain"]:
            prev_mc_std = [get_substance_name(x) for x in prev_conditions["must_contain"]]
            names_list = ", ".join([f"'{x}'" for x in prev_mc_std])
            query += f"AND s{step_index}.name IN [{names_list}] "
            query += f"AND {product_rel}.mainProduct = 'true' "
        # End bridge-compound constraint.
        

        for cpd in [ 'Nitrogen']:             # Exclude nitrogen as an intermediate compound.
            query += f"AND s{step_index}.name <> \"{cpd}\" "

        if "conversion_limitation" in cond and cond["conversion_limitation"] > 0:
            query += f"""
            AND {reactant_rel}.conversion IS NOT NULL
            AND toFloat(REPLACE({reactant_rel}.conversion, '%', '')) > {cond['conversion_limitation']}
            """
        if "yield_limitation" in cond and cond["yield_limitation"] > 0:
            query += f"""
            AND {product_rel}.yield IS NOT NULL
            AND toFloat(REPLACE({product_rel}.yield, '%', '')) > {cond['yield_limitation']}
            """
        if "selectivity_limitation" in cond and cond["selectivity_limitation"] > 0:
            query += f"""
            AND {product_rel}.selectivity IS NOT NULL
            AND toFloat(REPLACE({product_rel}.selectivity, '%', '')) > {cond['selectivity_limitation']}
            """
        if "score" in cond and cond["score"] > 0:
            query += f"""
            AND {product_rel}.score IS NOT NULL
            AND toFloat({product_rel}.score) > {cond['score']}
            """
        if cond.get("isHeterogeneous", "") == "true":
            query += f"AND {step_label}.isHeterogeneous = 'true' "
        if "include_ids" in cond and cond["include_ids"]:
            for id_val in cond['include_ids']:
                query += f"AND {step_label}.ID = '{id_val}' "
        if "exclude_ids" in cond and cond["exclude_ids"]:
            for id_val in cond["exclude_ids"]:
                query += f"AND {step_label}.ID <> '{id_val}' "
        if "must_contain" in cond:
            for cpd in cond["must_contain"]:
                query += f"AND '{get_substance_name(cpd)}' IN {jlist(f'{step_label}.reactantList')} "

        if "must_contain_products" in cond:
            for cpd in cond["must_contain_products"]:
                query += f"AND '{get_substance_name(cpd)}' IN {jlist(f'{step_label}.productList')} "

        if "must_not_contain" in cond:
            for cpd in cond["must_not_contain"]:
                std_name = get_substance_name(cpd)
                query += f"""
                AND NOT '{std_name}' IN {jlist(f'{step_label}.reactantList')}
                AND NOT '{std_name}' IN {jlist(f'{step_label}.solventList')}
                AND ALL(g IN apoc.convert.fromJsonList({step_label}.gases) WHERE g.name IS NULL OR g.name <> '{std_name}')
                """

        if "exclude_products" in cond:
            for cpd in cond["exclude_products"]:
                query += f"AND NOT '{get_substance_name(cpd)}' IN {jlist(f'{step_label}.productList')} "

                
        if "max_product_names" in cond and cond['max_product_names'] not in [False, 0]:
            query += f"AND size({jlist(f'{step_label}.productList')}) < {cond['max_product_names']+1} "
        if "max_reactant_names" in cond and cond['max_reactant_names'] not in [False, 0]:
            query += f"AND size({jlist(f'{step_label}.reactantList')}) < {cond['max_reactant_names']+1} "
        if cond.get("isSingleReaction", "") == "true":
            query += f"AND {step_label}.singleReaction = 'true' "
        if "temperature" in cond and cond["temperature"] != False:
            tmin, tmax = cond["temperature"]
            query += f"""
            AND toFloat(REPLACE({step_label}.temperature, ' K', '')) >= {tmin}
            AND toFloat(REPLACE({step_label}.temperature, ' K', '')) <= {tmax}
            """
        if "citations" in cond and cond["citations"] > 0:
            query += f"AND toInteger({step_label}.citations) >= {cond['citations']} "
        if "impactEffect" in cond and cond["impactEffect"] > 0:
            query += f"AND toFloat({step_label}.impactEffect) >= {cond['impactEffect']} "
        if "reactionType" in cond and cond["reactionType"]:
            query += f"AND {step_label}.reactionType = '{cond['reactionType']}' "
        if "catalyst" in cond and cond["catalyst"]:
            query += f"AND {step_label}.catalystName CONTAINS '{cond['catalyst']}' "
        if "pubDate" in cond and cond["pubDate"]:
            query += f"AND toInteger(substring({step_label}.pubDate, 0, 4)) >= {cond['pubDate']} "

    len_steps = len(step_conditions)
    reaction_expr = ", ".join([f"apoc.map.removeKeys(m{i}, ['new_pubDate'])" for i in range(len_steps)]).replace("m0", "m")
    # Final return clause
    if len(step_conditions) == 1:
        return_clause = f"""
        WITH DISTINCT [r.score] AS score, 
                      [{reaction_expr}] AS reaction, 
                      [j.journalName] AS journal, 
                      [s.name, {jlist('m.reactantList')}] AS path, 
                      [r.selectivity, r.yield] AS relation
        RETURN score, reaction, journal, path, relation
        LIMIT {path_limitation}
        """
    elif len(step_conditions) == 2:
        return_clause = f"""
        WITH DISTINCT [r.score, r2.score] AS score, 
                      [{reaction_expr}] AS reaction, 
                      [j.journalName, j1.journalName] AS journal, 
                      [s1.name, {jlist('m1.reactantList')}] AS path, 
                      [r2.selectivity, r1.conversion, r.selectivity, r.yield] AS relation
        RETURN score, reaction, journal, path, relation
        LIMIT {path_limitation}
        """
    elif len(step_conditions) == 3:
        return_clause = f"""
        MATCH (m2)<-[r5:is_reactant_of]-(s3:Compound)
        WITH DISTINCT [r.score, r2.score, r4.score] AS score, 
                      [{reaction_expr}] AS reaction, 
                      [j.journalName, j1.journalName, j2.journalName] AS journal, 
                      [s1.name, s2.name, {jlist('m2.reactantList')}] AS path,  
                      [r4.selectivity, r3.conversion, r2.selectivity, r1.conversion, r.selectivity, r.yield] AS relation
        RETURN score, reaction, journal, path, relation
        LIMIT {path_limitation}
        """
    elif len(step_conditions) == 4:
        return_clause = f"""
        MATCH (m3)<-[r7:is_reactant_of]-(s4:Compound)
        WITH DISTINCT [r.score, r2.score, r4.score, r6.score] AS score, 
                      [{reaction_expr}] AS reaction, 
                      [j.journalName, j1.journalName, j2.journalName, j3.journalName] AS journal, 
                      [s1.name, s2.name, s3.name, {jlist('m3.reactantList')}] AS path,  
                      [r6.selectivity, r5.conversion, r4.selectivity, r3.conversion, r2.selectivity, r1.conversion, r.selectivity, r.yield] AS relation
        RETURN score, reaction, journal, path, relation
        LIMIT {path_limitation}
        """
    else:
        return_clause = "WITH DISTINCT * RETURN *"


    reactant_sets = [f"{jlist('m.reactantList')}"]
    if len(step_conditions) >= 2:
        reactant_sets.append(f"{jlist('m1.reactantList')}")
    if len(step_conditions) >= 3:
        reactant_sets.append(f"{jlist('m2.reactantList')}")
    if len(step_conditions) >= 4:
        reactant_sets.append(f"{jlist('m3.reactantList')}")


    query += return_clause
    return query

def apply_all_conditions(step_conditions, all_conditions):
    'Merge shared flags, bibliographic constraints, scores, and product-count limits into step conditions. Existing nonempty heterogeneous and single-reaction flags take precedence. Percentage thresholds and temperatures are not inherited by this implementation.'
    updated_conditions = []
    if not all_conditions:
        # Return step conditions when no shared conditions are supplied.
        return step_conditions

    # Use the first shared-condition dictionary.
    ac = all_conditions[0]
    if "max_product_names" not in ac:
        ac["max_product_names"] = 8
        
    for step in step_conditions:
        # Copy the current step conditions.
        new_step = step.copy()


        # isHeterogeneous
        # Use shared values for missing or empty flags.
        if "isHeterogeneous" not in new_step or not new_step["isHeterogeneous"]:
            new_step["isHeterogeneous"] = ac.get("isHeterogeneous", "")

        # isSingleReaction
        if "isSingleReaction" not in new_step or not new_step["isSingleReaction"]:
            new_step["isSingleReaction"] = ac.get("isSingleReaction", "")

        # Copy shared bibliographic, score, and product-count constraints.
        for k in ["citations","impactEffect","pubDate","score", "max_product_names"]:
            if k in ac:
                new_step[k] = ac[k]

        updated_conditions.append(new_step)

    return updated_conditions

def _normalized_names(names: List[str]) -> List[str]:
    normalized = []
    for name in names or []:
        std_name = get_substance_name(name)
        if std_name:
            normalized.append(std_name)
    return normalized

def _normalized_name_set(names: List[str]) -> set[str]:
    return set(_normalized_names(names))

def generate_split_query(
    target_name: str, 
    n_steps: int, 
    group_pattern: List[str],      # e.g. ['H', 'T', 'H']
    current_cutoff: float,         
    step_conditions_data: List[dict], 
    path_limitation: int,          
    search_config: Dict[str, Any], 
    main_product: str = "true"
) -> Optional[str]:
    'Build a Cypher query with per-step branch-and-bound pruning.'

    # Read branch-and-bound settings.
    max_potential_map = search_config["max_potential_map"]
    score_h_min = search_config["score_h_min"]
    
    # Skip groups whose upper bound is below the cutoff.
    current_group_potential = sum(max_potential_map.get(p, 0) for p in group_pattern)
    if current_group_potential < current_cutoff:
        return None

    # Cypher list helper.
    def jlist(field): 
        return f"coalesce(apoc.convert.fromJsonList({field}), [])"

    # Build the query.
    return _build_query_core(
        target=target_name,
        n=n_steps,
        pattern=group_pattern,
        cutoff=current_cutoff,
        max_potential_map=max_potential_map,
        h_min=score_h_min,
        conditions=step_conditions_data,
        config=search_config,
        main_prod=main_product,
        limit=path_limitation,
        jlist=jlist,
    )

def _build_query_core(target, n, pattern, cutoff, max_potential_map, h_min, conditions, config, main_prod, limit, jlist):
    lines = []
    
    # Read per-query settings.
    banned_intermediates = _normalized_names(config.get("banned_intermediates", []))
    mutex_groups = config.get("mutex_groups", [])
    avoid_loops = config.get("avoid_loops", True)
    exclude_target = config.get("exclude_target_in_reactants", True)
    global_temp = config.get("global_temperature_range", [None, None])

    # Build MATCH clauses and pruning conditions.
    lines.append(f"MATCH (s:Compound {{name: '{get_substance_name(target)}'}})")
    
    m_vars = []   
    s_vars = ["s"] 
    
    # Retain variables that must remain in scope across WITH clauses.
    collected_vars = ["s"] 
    
    for i in range(n):
        m_curr = "m" if i == 0 else f"m{i}"
        step_type = pattern[i]
        cond = conditions[i] if i < len(conditions) else {}
        impact_expr = f"coalesce(toFloat({m_curr}.impactEffect), 0.0)"
        
        # Relationship variables for this step.
        rel_prod = "r" if i == 0 else f"r{i*2}"
        rel_react = "" if i == 0 else f"r{i*2-1}"
        
        # Match the current step.
        if i == 0:
            lines.append(f"MATCH (s)<-[{rel_prod}:has_product]-({m_curr}:Reaction)")
            # Track variables in scope.
            collected_vars.append(rel_prod)
            collected_vars.append(m_curr)
        else:
            m_prev = "m" if i == 1 else f"m{i-1}"
            s_prev = "s" if i == 0 else f"s{i}" 
            lines.append(f"MATCH ({m_prev})<-[{rel_react}:is_reactant_of]-({s_prev}:Compound)<-[{rel_prod}:has_product]-({m_curr}:Reaction)")
            s_vars.append(s_prev)
            # Track newly introduced variables in order.
            collected_vars.append(rel_react)
            collected_vars.append(s_prev)
            collected_vars.append(rel_prod)
            collected_vars.append(m_curr)

        m_vars.append(m_curr)
        
        # Start a WHERE clause after MATCH.
        prefix = "WHERE"
        
        # Apply single-step filters.
        
        # Exclude specified intermediate compounds.
        if i > 0 and banned_intermediates:
            prev_required_intermediates = _normalized_name_set(conditions[i - 1].get("must_contain", []))
            effective_banned_intermediates = [
                name for name in banned_intermediates
                if name not in prev_required_intermediates
            ]
            check_str = str(effective_banned_intermediates)
        else:
            effective_banned_intermediates = []

        if i > 0 and effective_banned_intermediates:
            lines.append(f"{prefix} NOT {s_vars[-1]}.name IN {check_str}")
            prefix = "AND"

        # Require the specified bridge compound as the product of this upstream step.
        if i > 0:
            prev_must_contain = _normalized_names(conditions[i - 1].get("must_contain", []))
            if prev_must_contain:
                lines.append(f"{prefix} {s_vars[-1]}.name IN {prev_must_contain}")
                lines.append(f"AND {rel_prod}.mainProduct = 'true'")
                prefix = "AND"

        # Apply H/T score constraints.
        if step_type == 'H':
            lines.append(f"{prefix} {impact_expr} > {h_min}")
        else:
            lines.append(f"{prefix} {impact_expr} <= {h_min}")
        prefix = "AND"

        # Apply required and excluded reactant constraints.
        must_contain = _normalized_names(cond.get("must_contain", []))
        for compound_name in must_contain:
            lines.append(f"{prefix} '{compound_name}' IN {jlist(f'{m_curr}.reactantList')}")
            prefix = "AND"

        must_contain_products = _normalized_names(cond.get("must_contain_products", []))
        for compound_name in must_contain_products:
            lines.append(f"{prefix} '{compound_name}' IN {jlist(f'{m_curr}.productList')}")
            prefix = "AND"

        must_not_contain = _normalized_names(cond.get("must_not_contain", []))
        for compound_name in must_not_contain:
            lines.append(f"{prefix} NOT '{compound_name}' IN {jlist(f'{m_curr}.reactantList')}")
            lines.append(f"AND NOT '{compound_name}' IN {jlist(f'{m_curr}.solventList')}")
            lines.append(
                f"AND ALL(g IN {jlist(f'{m_curr}.gases')} "
                f"WHERE g.name IS NULL OR g.name <> '{compound_name}')"
            )
            prefix = "AND"

        exclude_products = _normalized_names(cond.get("exclude_products", []))
        for compound_name in exclude_products:
            lines.append(f"{prefix} NOT '{compound_name}' IN {jlist(f'{m_curr}.productList')}")
            prefix = "AND"

        include_ids = cond.get("include_ids", [])
        for reaction_id in include_ids:
            lines.append(f"{prefix} {m_curr}.ID = '{reaction_id}'")
            prefix = "AND"

        exclude_ids = cond.get("exclude_ids", [])
        for reaction_id in exclude_ids:
            lines.append(f"{prefix} {m_curr}.ID <> '{reaction_id}'")
            prefix = "AND"

        # Calculate the upper bound for the remaining steps.
        future_potential = sum(max_potential_map.get(pattern[k], 0) for k in range(i + 1, n))

        # Apply the initial-step bound.
        if i == 0:
            min_val = cutoff - future_potential
            lines.append(f"AND {impact_expr} > {max(0, min_val)}")

        # Apply temperature and property filters.
        step_temp = cond.get("temperature")
        t_min = step_temp[0] if step_temp else global_temp[0]
        t_max = step_temp[1] if step_temp else global_temp[1]
        
        if t_min is not None:
            lines.append(f"AND toFloat(REPLACE({m_curr}.temperature, ' K', '')) >= {t_min}")
        if t_max is not None:
            lines.append(f"AND toFloat(REPLACE({m_curr}.temperature, ' K', '')) <= {t_max}")

        if cond.get("isHeterogeneous", "") == "true":
            lines.append(f"AND {m_curr}.isHeterogeneous = 'true'")

        if cond.get("isSingleReaction", "") == "true":
            lines.append(f"AND {m_curr}.singleReaction = 'true'")

        if cond.get("citations", 0) > 0:
            lines.append(f"AND toInteger({m_curr}.citations) >= {cond['citations']}")

        if cond.get("impactEffect", 0) > 0:
            lines.append(f"AND {impact_expr} >= {cond['impactEffect']}")

        if cond.get("reactionType"):
            lines.append(f"AND {m_curr}.reactionType = '{cond['reactionType']}'")

        if cond.get("catalyst"):
            lines.append(f"AND {m_curr}.catalystName CONTAINS '{cond['catalyst']}'")

        if cond.get("pubDate"):
            lines.append(f"AND toInteger(substring({m_curr}.pubDate, 0, 4)) >= {cond['pubDate']}")

        if cond.get("max_product_names") not in [None, False, 0]:
            lines.append(f"AND size({jlist(f'{m_curr}.productList')}) < {cond['max_product_names'] + 1}")

        if cond.get("max_reactant_names") not in [None, False, 0]:
            lines.append(f"AND size({jlist(f'{m_curr}.reactantList')}) < {cond['max_reactant_names'] + 1}")
        
        if main_prod == "true":
            lines.append(f"AND {rel_prod}.mainProduct = 'true'")

        # Prune after matching and filtering each step.
        # Check the accumulated score against the remaining upper bound.
        # Condition: current prefix sum + future potential >= cutoff.
        if i < n - 1:
            # Sum the scores of the steps matched so far.
            current_impact_vars = [f"coalesce(toFloat({m}.impactEffect), 0.0)" for m in m_vars]
            sum_expr = " + ".join(current_impact_vars)

            remaining_steps = n - i - 1
            future_bound_field = f"ub_if_{remaining_steps}"

            # Use the most promising upstream reactant to bound the remaining score.
            lines.append(f"OPTIONAL MATCH ({m_curr})<-[:is_reactant_of]-(future_c:Compound)")

            # Build the WITH clause.
            vars_str = ", ".join(collected_vars)
            lines.append(
                f"WITH {vars_str}, "
                f"coalesce(max(coalesce(future_c.{future_bound_field}, 0.0)), 0.0) AS future_bound"
            )
            
            # Apply the pruning condition.
            lines.append(f"WHERE ({sum_expr} + future_bound >= {cutoff})")

    # Apply whole-path constraints.
    
    has_post_with_constraints = False

    # Check mutually exclusive groups.
    if mutex_groups:
        parts = []
        for m in m_vars:
            parts.append(f"{jlist(f'{m}.reactantList')}")
            parts.append(f"{jlist(f'{m}.solventList')}")
            parts.append(f"[g IN {jlist(f'{m}.gases')} | g.name]")
        
        union_str = " + ".join(parts)
        # Preserve the variables collected by the MATCH clauses.
        # WITH * also retains variables introduced by APOC.
        lines.append(f"WITH *, apoc.coll.toSet(apoc.coll.flatten([{union_str}])) AS all_chemicals")
        
        mutex_prefix = "WHERE"
        for group in mutex_groups:
            subs = group.get("substances", [])
            if len(subs) >= 2:
                check_list = str(subs)
                lines.append(f"{mutex_prefix} NOT ALL(x IN {check_list} WHERE x IN all_chemicals)")
                mutex_prefix = "AND"
                has_post_with_constraints = True

    # Remove repeated compounds along the path.
    if avoid_loops:
        s_names = [f"{s}.name" for s in s_vars]
        s_names_str = ", ".join(s_names)
        
        # Choose the clause prefix.
        if mutex_groups:
            loop_prefix = "AND" if has_post_with_constraints else "WHERE"
        else:
            # Without mutually exclusive groups, continue the existing WHERE clause.
            # The last step does not add a pruning WITH clause.
            # Append the remaining constraints with AND.
            loop_prefix = "AND"
            
        lines.append(f"{loop_prefix} size(apoc.coll.toSet([{s_names_str}])) = {len(s_vars)}")

    # Exclude the target product from starting materials.
    if exclude_target:
        for m in m_vars:
            lines.append(f"AND NOT '{get_substance_name(target)}' IN {jlist(f'{m}.reactantList')}")
            
    # Check the final total against the cutoff.
    # Recheck the total after the final step.
    all_impacts = [f"coalesce(toFloat({m}.impactEffect), 0.0)" for m in m_vars]
    lines.append(f"AND ({' + '.join(all_impacts)} >= {cutoff})")

    # Build the RETURN clause.
    return _generate_final_return(lines, n, limit, jlist)

def _generate_final_return(lines, n, limit, jlist):
    'Build the publication joins and RETURN clause.'
    for k in range(n):
        m_label = "m" if k == 0 else f"m{k}"
        j_label = "j" if k == 0 else f"j{k}"
        lines.append(f"MATCH ({m_label})<-[:is_article_of]-(:Article)-[:is_published_in]->({j_label}:Journal)")

    m_list = ["m"] + [f"m{k}" for k in range(1, n)]
    impact_vars = [f"coalesce(toFloat({x}.impactEffect), 0.0)" for x in m_list]
    sort_expression = " + ".join(impact_vars)

    with_clause = f"WITH *, {sort_expression} AS total_sort_key"

    if n == 1:
        return_clause = f"""
        {with_clause}
        RETURN total_sort_key, 
               [m.impactEffect] AS score, 
               [m] AS reaction, 
               [j.journalName] AS journal, 
               [s.name, {jlist('m.reactantList')}] AS path, 
               [r.selectivity, r.yield] AS relation
        ORDER BY total_sort_key DESC
        LIMIT {limit}
        """
    elif n == 2:
        return_clause = f"""
        {with_clause}
        RETURN total_sort_key,
               [m.impactEffect, m1.impactEffect] AS score, 
               [m, m1] AS reaction, 
               [j.journalName, j1.journalName] AS journal, 
               [s1.name, {jlist('m1.reactantList')}] AS path, 
               [r2.selectivity, r1.conversion, r.selectivity, r.yield] AS relation
        ORDER BY total_sort_key DESC
        LIMIT {limit}
        """
    elif n == 3:
        return_clause = f"""
        {with_clause}
        RETURN total_sort_key,
               [m.impactEffect, m1.impactEffect, m2.impactEffect] AS score, 
               [m, m1, m2] AS reaction, 
               [j.journalName, j1.journalName, j2.journalName] AS journal, 
               [s1.name, s2.name, {jlist('m2.reactantList')}] AS path,  
               [r4.selectivity, r3.conversion, r2.selectivity, r1.conversion, r.selectivity, r.yield] AS relation
        ORDER BY total_sort_key DESC
        LIMIT {limit}
        """
    elif n == 4:
        return_clause = f"""
        {with_clause}
        RETURN total_sort_key,
               [m.impactEffect, m1.impactEffect, m2.impactEffect, m3.impactEffect] AS score, 
               [m, m1, m2, m3] AS reaction, 
               [j.journalName, j1.journalName, j2.journalName, j3.journalName] AS journal, 
               [s1.name, s2.name, s3.name, {jlist('m3.reactantList')}] AS path,  
               [r6.selectivity, r5.conversion, r4.selectivity, r3.conversion, r2.selectivity, r1.conversion, r.selectivity, r.yield] AS relation
        ORDER BY total_sort_key DESC
        LIMIT {limit}
        """
    else:
        return_clause = f"{with_clause} RETURN * LIMIT {limit}"

    lines.append(return_clause)
    return "\n".join(lines)

REACTION_QUERY = """
UNWIND $reaction_ids AS rid
MATCH (r:Reaction {ID:rid})
CALL {
    WITH r
    OPTIONAL MATCH (a:Article)-[:is_article_of]->(r)
    RETURN collect(properties(a)) AS articles
}
CALL {
    WITH r
    OPTIONAL MATCH (cat:Catalyst)-[:is_catalyst_of]->(r)
    RETURN collect(properties(cat)) AS catalysts
}
CALL {
    WITH r
    OPTIONAL MATCH (c:Compound)-[e:is_reactant_of]->(r)
    RETURN collect({name:c.name, props:properties(e)}) AS reactants
}
CALL {
    WITH r
    OPTIONAL MATCH (r)-[e:has_product]->(c:Compound)
    RETURN collect({name:c.name, props:properties(e)}) AS products
}
RETURN properties(r) AS node, articles, catalysts, reactants, products
ORDER BY rid
"""

PATH_QUERY = """
UNWIND $reaction_ids AS rid
MATCH (photo:Reaction {ID:rid, reactionType:'photocatalysis'})
MATCH (feed:Compound {name:$feed})-[:is_reactant_of]->(photo)
MATCH (h:Compound {name:$hydrogen})-[:is_reactant_of]->(photo)
MATCH (photo)-[p:has_product]->(co:Compound {name:$intermediate})
      -[:is_reactant_of]->(thermal:Reaction {reactionType:'thermal catalysis'})
      -[t:has_product]->(methanol:Compound {name:$target})
MATCH (h)-[:is_reactant_of]->(thermal)
WHERE thermal.ID IN $reaction_ids
  AND toLower(toString(p.mainProduct)) = 'true'
  AND toLower(toString(t.mainProduct)) = 'true'
  AND COUNT { (:Compound)-[:is_reactant_of]->(photo) } = 2
  AND COUNT { (:Compound)-[:is_reactant_of]->(thermal) } = 2
  AND COUNT { (photo)-[:has_product]->(:Compound) } <= 8
  AND COUNT { (thermal)-[:has_product]->(:Compound) } <= 8
  AND toFloat(replace(photo.temperature, ' K', '')) IS NOT NULL
  AND toFloat(replace(thermal.temperature, ' K', '')) IS NOT NULL
RETURN DISTINCT photo.ID AS photo_ID, thermal.ID AS thermal_ID
ORDER BY photo_ID, thermal_ID
"""
