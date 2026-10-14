# CatKG

CatKG organizes catalytic reaction information from published literature in a knowledge graph to support the search and comparison of multistep pathways. Our earlier study, [*Synergizing a knowledge graph and large language model for relay catalysis pathway recommendation*](https://doi.org/10.1093/nsr/nwaf271) (Fu et al., *National Science Review*, 12(8), nwaf271, 2025), focused on thermocatalysis and combined LLM-assisted literature extraction, knowledge-graph construction, chemistry-based pathway scoring, and LLM-assisted pathway analysis. This approach connects reactions through shared intermediates and allows the reaction evidence to be traced to the original publications.

MultiCat-KG extends this approach to photocatalysis, electrocatalysis, and thermocatalysis. It represents reactants, products, catalysts, reaction conditions, performance, and literature sources in a common graph, including illumination and electrochemical parameters where relevant. The graph supports the retrieval and ranking of candidate tandem pathways across these catalytic modes.

This repository provides the extraction prompts, data-processing and graph-construction code, pathway queries and scoring functions, and a small sample dataset for MultiCat-KG. It is intended as a practical reference for researchers interested in extracting catalytic data from literature, cleaning reaction records, and exploring candidate pathways. The supplied workflow shows how article information is converted into reaction records, connected in a graph, and used for pathway ranking.

Repository: [chenggroup/catkg](https://github.com/chenggroup/catkg).

## Literature extraction

This component converts catalytic information in article text, captions, and available tables into structured reaction records. It extracts reactants, products, catalysts, conditions, and performance, and associates them with the corresponding reported reaction.

`extraction/` contains `extract_reactions.py`, the Gemini 2.5 Flash configuration in `model_settings.json`, and the extraction prompts in `prompts.json`. The script performs coarse and detailed extraction, each followed by a review step.

### Extract a new article

Set `GOOGLE_API_KEY` to a Google API key with model access, then run:

```shell
python -m extraction.extract_reactions --input data/fulltexts/photocatalysis_01.json --output outputs_live/photocatalysis_01_extraction.json --live
```

The supplied workflow uses saved extraction results by default; a new Gemini extraction requires internet access and API credentials.

## Data processing and graph construction

Reaction descriptions, compound names, and units vary across articles. This component normalizes the extracted records and organizes them in a knowledge graph. Consistent compound names allow reactions to connect through common intermediates, while links to articles, journals, and publishers retain their literature sources.

`data_processing/` contains the normalization and checking functions. `clean_reactions.py` runs the cleaning steps, including checks on reaction completeness and performance values; `build_graph.py` constructs reaction, compound, catalyst, article, journal, and publisher nodes and their relationships.

The default workflow saves the node and relationship tables to `outputs/graph/`. There are two ways to use the sample graph in Neo4j:

- **Restore the supplied dump:** Load `data/neo4j/neo4j.dump` into a stopped, empty Neo4j instance using `neo4j-admin database load`, then start the instance. The dump contains the prepared sample graph, so rebuilding it from the reaction records is unnecessary.
- **Build from the cleaned records:** Run the workflow to generate `outputs/graph/`, start a dedicated empty database, set the connection variables, and run `python data_processing/neo4j_graph.py --import-demo` from the repository directory. The script imports the nodes and relationships and checks the retrieved path scores against the offline results.

Both approaches support the Neo4j retrieval and scoring example. See the [Neo4j setup, dump restoration, and import instructions](data_processing/README.md#neo4j-demonstration). The supplied dump uses Neo4j Community 5.20.0 and Java 17, installed separately from the Python dependencies.

## Path recommendation

This component searches for reaction sequences in which a product of one step can serve as a reactant in the next. It ranks candidate paths using reported reaction performance, selected condition matches, and literature metadata, allowing researchers to compare candidates and inspect the supporting reaction records.

`path_recommendation/` contains `queries.py` for graph queries, `recommendation.py` for path retrieval and ranking, and `scoring.py` for single-reaction scores, adjacent-step compatibility, and whole-path contributions.

`scoring_and_ranking_example.py` prints each score contribution for the highest-ranked sample path, **284292 → 311263**, whose total score is **3702.93**. It supports saved records and read-only Neo4j retrieval. `scoring_and_ranking_example.csv` provides the first 20 ranked paths for comparison. See the [scoring rules](path_recommendation/README.md).

## Sample data

The sample data provide a compact example of the workflow from article inputs and saved extraction results to cleaned reactions, graph construction, and pathway ranking. They allow readers to run the code and inspect the reaction information alongside its source articles.

`data/` contains 20 photocatalytic CO2 + H2-to-CO reactions and 20 thermocatalytic CO + H2-to-methanol reactions, together with their 40 source articles. The reaction records come from the cleaned **471,436-reaction Neo4j graph used in the paper**. The workflow reconstructs the sample graph and recalculates scores for its 400 connected two-step paths.

| File or folder | Contents and purpose |
|---|---|
| `fulltexts/` | Forty article JSON files used as extraction inputs |
| `samples.csv` | Correspondence between reaction IDs, DOIs, and article files |
| `extracted_reactions.json` | Saved extraction results for running the workflow |
| `compound_names.json` | Compound aliases and common names used to connect reactions |
| `graph_metadata.json` | Article, journal, publisher, and reaction metadata used to build the graph |
| `neo4j/neo4j.dump`, `neo4j/graph_info.json` | Sample graph dump, Neo4j version, checksum, and graph counts |

## Run the workflow

Use **Python 3.10 or later**; the workflow has been checked with Python 3.10.13. Python is installed separately; `requirements.txt` lists the Python packages. An isolated environment is optional.

From the repository directory, install the dependencies and start JupyterLab:

```shell
python -m pip install -r requirements.txt
python -m jupyterlab
```

Open `literature_to_tandem_paths.ipynb`, select the environment containing these dependencies, and run the cells. The notebook follows the four parts above. New Gemini extraction and Neo4j retrieval are optional cells.

To process the saved reaction records, construct the graph, and rank paths from the command line:

```shell
python process_reactions_and_rank_paths.py
python path_recommendation/scoring_and_ranking_example.py
```

Cleaned reactions, graph tables, rankings, and score details are written to `outputs/`. To calculate scores from Neo4j, use `python path_recommendation/scoring_and_ranking_example.py --source neo4j` after following the [connection instructions](data_processing/README.md#neo4j-demonstration).
