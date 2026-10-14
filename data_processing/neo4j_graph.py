"""Import/retrieve the small demo; refuse writes to an unrelated database."""
import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from neo4j import GraphDatabase, READ_ACCESS
from data_processing.clean_reactions import read, write
from path_recommendation.recommendation import retrieve_and_rank, sample_reaction_ids, compound_names, formatted_rankings

NODE_TABLES = [('Article','article_nodes','article_ID'),('Reaction','reaction_nodes','ID'),
               ('Catalyst','catalyst_nodes','demoCatalystID'),('Compound','compound_nodes','name'),
               ('Journal','journal_nodes','journalName'),('Publisher','publisher_nodes','publisher')]
SCOPE = 'multicatkg-40-demo'


def props(row):
    return {k:json.dumps(v,ensure_ascii=False) if isinstance(v,dict) or isinstance(v,list) and any(isinstance(x,(dict,list)) for x in v) else v
            for k,v in row.items() if v is not None}


def import_graph(session, graph):
    existing = session.run('MATCH (n) RETURN count(n) AS n').single()['n']
    expected = {x['ID'] for x in graph['reaction_nodes']}
    if len(expected) != 40:raise RuntimeError('Import requires exactly 40 demo reactions')
    if existing:
        current = session.run('MATCH (r:Reaction) RETURN r.ID AS id,r.demoScope AS scope').data()
        if {x['id'] for x in current} != expected or any(x['scope'] != SCOPE for x in current):
            raise RuntimeError('Refusing import: target is not an empty database or this dedicated demo')
    for label,table,key in NODE_TABLES:
        session.run(f'CREATE CONSTRAINT demo_{label.lower()} IF NOT EXISTS FOR (n:{label}) REQUIRE n.{key} IS UNIQUE').consume()
        rows = []
        for r in graph[table]:
            r = dict(r)
            if label=='Catalyst':r[key]=str(r['catalystName'])+'||'+str(r['reactionID'])
            if label=='Reaction':r['demoScope']=SCOPE
            rows.append({'key':r[key],'props':props(r)})
        session.run(f'UNWIND $rows AS row MERGE (n:{label} {{{key}:row.key}}) SET n += row.props',rows=rows).consume()
    for role,rel in [('reactant','is_reactant_of'),('solvent','is_solvent_of'),('product','has_product')]:
        rows=[{'rid':e['r'],'name':e['c'],'props':props(e['p'])} for e in graph['compoundToReaction_relation'] if e['t']==role]
        match='(r)-[e:has_product]->(c)' if role=='product' else f'(c)-[e:{rel}]->(r)'
        session.run(f'UNWIND $rows AS row MATCH (r:Reaction {{ID:row.rid}}),(c:Compound {{name:row.name}}) MERGE {match} SET e += row.props',rows=rows).consume()
    session.run('UNWIND $rows AS row MATCH (a:Article {article_ID:row.aid}),(r:Reaction {ID:row.rid}) MERGE (a)-[:is_article_of]->(r)',
                rows=[{'aid':str(e['aid']),'rid':str(e['rid'])} for e in graph['reactionToArticle_relation']]).consume()
    session.run('MATCH (c:Catalyst),(r:Reaction) WHERE c.reactionID=r.ID MERGE (c)-[:is_catalyst_of]->(r)').consume()
    session.run('UNWIND $rows AS row MATCH (a:Article {article_ID:row.aid}),(j:Journal {journalName:row.journalName}) MERGE (a)-[:is_published_in]->(j)',rows=graph['articleToJournal_relation']).consume()
    session.run('UNWIND $rows AS row MATCH (j:Journal {journalName:row.journalName}),(p:Publisher {publisher:row.publisher}) MERGE (j)-[:is_published_by]->(p)',rows=graph['journalToPublisher_relation']).consume()


def query_graph(session):
    scope=session.run('MATCH (r:Reaction) RETURN count(r) AS n,sum(CASE WHEN r.demoScope=$scope THEN 1 ELSE 0 END) AS marked',scope=SCOPE).single()
    if scope['n']!=40 or scope['marked']!=40:
        raise RuntimeError('Query target must be the dedicated 40-reaction demonstration database')
    ranked, records, stats = session.execute_read(retrieve_and_rank, sample_reaction_ids(ROOT), compound_names(ROOT))
    paths = formatted_rankings(ranked, records)
    counts=session.run('MATCH (n) RETURN labels(n)[0] AS label,count(n) AS count ORDER BY label').data()
    relations=session.run('MATCH ()-[e]->() RETURN type(e) AS relation,count(e) AS count ORDER BY relation').data()
    expected_links={'is_article_of':40,'is_published_in':40}
    observed={x['relation']:x['count'] for x in relations}
    if any(observed.get(k)!=v for k,v in expected_links.items()):raise RuntimeError('Incomplete article provenance links in demo graph')
    return paths,{'source':'Neo4j queried demo graph','nodes':counts,'relations':relations,'ranking':stats,'retrieved_path_pairs':stats['retrieved_pairs']}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--uri',default=os.getenv('NEO4J_URI','bolt://localhost:7687'))
    ap.add_argument('--database',default=os.getenv('NEO4J_DATABASE','neo4j'))
    ap.add_argument('--no-auth', action='store_true', help='Use only if database authentication is disabled')
    ap.add_argument('--import-demo',action='store_true')
    ap.add_argument('--graph',type=Path,default=ROOT/'outputs/graph')
    ap.add_argument('--output',type=Path,default=ROOT/'outputs/neo4j_validation')
    args=ap.parse_args()
    if not (ROOT / 'outputs/path_rankings.json').exists():
        ap.error('Run python process_reactions_and_rank_paths.py first to generate the offline reference ranking')
    password=os.getenv('NEO4J_PASSWORD')
    if not args.no_auth and not password:
        import getpass
        password = getpass.getpass('Neo4j password: ')
    auth = None if args.no_auth else (os.getenv('NEO4J_USERNAME', 'neo4j'), password)
    with GraphDatabase.driver(args.uri,auth=auth) as driver:
        if args.import_demo:
            graph={p.stem:read(p) for p in args.graph.glob('*.json')}
            with driver.session(database=args.database) as session:import_graph(session,graph)
        with driver.session(database=args.database,default_access_mode=READ_ACCESS) as session:
            paths,summary=query_graph(session)
    reference=read(ROOT/'outputs/path_rankings.json')
    def scores(p):return {(x['photo_ID'],x['thermal_ID']):x['total'] for x in p}
    if scores(paths)!=scores(reference):raise RuntimeError('Neo4j and offline path scores differ')
    write(args.output/'path_rankings.json',paths)
    summary['offline_score_match']=True
    write(args.output/'verification.json',summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
