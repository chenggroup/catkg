# Path recommendation

| File | Purpose |
|---|---|
| `queries.py` | Compound-name lookup, general backward and branch-and-bound queries, and parameterized sample queries |
| `scoring.py` | Single-reaction A1–A8, adjacent-step B1–B8, whole-path C1–C2, and ranking |
| `recommendation.py` | Retrieve connected sample paths offline or from Neo4j and recalculate scores |
| `scoring_and_ranking_example.py` | Print the calculation for the highest-ranked sample path |
| `scoring_and_ranking_example.csv` | Reference ranking for the first 20 sample paths |

Run from the package directory:

```shell
python path_recommendation/scoring_and_ranking_example.py
python path_recommendation/scoring_and_ranking_example.py --source neo4j
```

The default uses the supplied extraction results after cleaning. The Neo4j option reads the reaction, compound, catalyst, and article properties and connected paths from the database. The sample index supplies reaction IDs only; stored product scores are checked, and scores are recalculated from the retrieved properties. Both modes call the same scoring functions and print A1–A8, B1–B8, and C1–C2. For the supplied sample, reaction **284292 → 311263** ranks first at **3702.93**. This is a calculation example of a candidate path.

The output includes all ranked pairs in `path_rankings.json` and a query summary in `verification.json`, under `outputs/scoring_example/` or `outputs/neo4j_scoring/`. Use `--output DIRECTORY` to choose another location. Neo4j setup and connection variables are described in [data processing](../data_processing/README.md#neo4j-demonstration).

The general queries support other target compounds and condition limits and require APOC. The fixed sample queries require no APOC.

## Scoring rules and missing values

The score orders candidate paths by reported reaction performance, publication metadata, field completeness, and selected condition matches. The example recalculates the A/B/C contributions using the supplied implementation.

## Category A: one reaction

| Term | Implemented rule |
|---|---|
| A1 | `round(IF / 50.5 * 200, 2)`; not capped when IF exceeds 50.5 |
| A2 | Linear date score from 2000-01-01 to **2024-07-30**, capped at 100 for later dates; earlier/unparseable dates: 0 |
| A3 | Citation count when 0–100; **200 when count >100**; unparseable/negative: 0 |
| A4 | 20 for each of 10 available catalyst fields; at most 200 |
| A5 | Temperature field available: 100; missing: 0 |
| A6 | Nonempty gas list and every gas has a name and pressure: 100; otherwise 0 |
| A7 | Reaction time available: 50; missing: 0 |
| A8 | Sum of performance contributions over **all products** |

A4 fields are loading, support, preparation method, promoters, metal oxidation states, BET area, TOF, shape, pore size, and pore volume. Availability checks mostly compare values with `None`, `unknown`, empty lists/strings. A4 counts available catalyst fields. Light-source information and Faradaic efficiency are stored in the graph; the A formula uses the terms listed above.

For maximum reactant conversion `x`, product selectivity `s`, and yield `y`, all represented as fractions:

```text
product contribution = 400 * max(s*x, y) + 100 * (2*s + x + y)
A8 = sum(product contribution for every product)
A = A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8
```

Missing conversion/yield/selectivity contribute zero to their numeric variables. Yield `unknown` does not imply zero observed product formation. A known conversion can still contribute through `100*x` even for a product whose yield and selectivity are both unknown. Because A8 sums over products, the shared conversion contribution repeats for each product. Every product relationship receives the **same reaction-level total A**. 

Run cleaning before scoring so that performance inputs use the percentage format expected by the scoring function.

## Category B: adjacent reactions

| Term | Implemented rule |
|---|---|
| B1 | Same known `isHeterogeneous` flag: 50; this does not compare catalytic modalities |
| B2 | Same known catalyst name, case-insensitive: 200 |
| B3 | Same known support, case-insensitive: 100 |
| B4 | Equal nonempty sets of promoter names: 100 |
| B5 | Numeric temperatures: `30000 / (abs(T1-T2) + 30)` |
| B6 | At least one common known gas name: 100 |
| B7 | Largest `200*(1-abs(p1-p2)/(p1+p2))` over shared gases with usable pressures |
| B8 | At least one identical solvent name: 100 |

B5 has no cross-modality exemption. Equal temperatures earn 1000. For 973 K and 298 K, the bonus is `30000/705 = 42.553191...`; missing temperatures earn 0. Large gaps reduce a positive bonus, rather than applying a negative penalty. B7 expects upstream-normalized pressure units: the path function itself strips some unit suffixes but does not convert them.

## Category C and total

```text
C1 = 500 if every reactant of the first forward step is in the simple-compound list
C2 = 10*(15-n) if n <= 15 else 0
total = round(sum(A_step_scores) + sum(B_adjacent_pair_scores) + C1 + C2, 2)
```

`n` counts literal names from reactants, gases, and solvents over the whole path. Product names are not directly added. Aliases such as `CO2` and `Carbon Dioxide` can both be counted. Retrieval requires nonempty reactant matches before scoring.

The path implementation drops a route when the counted literal-name set contains both **`H2` and `O2`**. The check uses the stored literal names.

## Retrieval order

Backward queries store the target-side reaction first. C1 reads the last reaction in that list, which is the first forward step. Ranked results display the photocatalytic step followed by the thermocatalytic step.

The demonstration scores all 400 pairs. `--temperature-shortlist N` optionally keeps the N pairs with the smallest temperature differences before applying A+B+C. The general query functions in `path_recommendation/queries.py` support reaction and condition constraints.

`generate_split_query()` in `queries.py` provides branch-and-bound query construction for the alternative search method; the sample example uses A+B+C scoring.
