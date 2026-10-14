"""Run the extraction stages with Gemini using the supplied prompts.
"""
import argparse
import csv
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def format_article_content(content):
    """Format article text, captions, and tables for extraction."""
    pick = lambda x: x[0] if isinstance(x, list) and x else x
    text = ""
    for field in ("title", "abstract"):
        if content.get(field):
            text += f"- {field}: {pick(content[field])}\n"
    for field in ("paragraphs", "figureCaptions", "schemeCaptions", "tables"):
        if field in content:
            text += f"- {field}:\n"
            for item in content[field]:
                text += f"  - {item}\n"
    return text


def extract_reactions(obj):
    """Collect reaction entries from the coarse extraction output."""
    found = []
    if isinstance(obj, dict):
        if all(k in obj for k in ("id", "reactants", "products", "catalyst")):
            return [obj]
        for value in obj.values():
            found.extend(extract_reactions(value))
    elif isinstance(obj, list):
        for value in obj:
            found.extend(extract_reactions(value))
    return found


def build_request(system, user_text, config):
    return {"systemInstruction": {"parts": [{"text": system}]},
            "contents": [{"role": "user", "parts": [{"text": user_text}]}],
            "generationConfig": config["generationConfig"], "safetySettings": config["safetySettings"]}


def generate(system, user_text, config, audit, stage):
    key = os.getenv("GOOGLE_API_KEY", "").strip()
    if not key:
        raise RuntimeError("Set GOOGLE_API_KEY before invoking live extraction")
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{config['model']}:generateContent"
    payload = build_request(system, user_text, config)
    request = urllib.request.Request(url, data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "x-goog-api-key": key}, method="POST")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                result = json.load(response)
            break
        except urllib.error.HTTPError as error:
            if error.code in (429, 500, 502, 503, 504) and attempt < 2:
                time.sleep(2 ** attempt)
                continue
            raise RuntimeError(f"Gemini returned HTTP {error.code}") from None
    candidates = result.get("candidates", [])
    parts = candidates[0].get("content", {}).get("parts", []) if candidates else []
    text = "".join(p.get("text", "") for p in parts if not p.get("thought", False))
    if not text:
        raise RuntimeError(f"No text response for {stage}; extraction remains incomplete")
    parsed = json.loads(text)
    audit.append({"stage": stage, "model": config["model"], "returned_model_version": result.get("modelVersion"),
        "system_prompt_sha256": hashlib.sha256(system.encode("utf-8")).hexdigest(),
        "user_message_sha256": hashlib.sha256(user_text.encode("utf-8")).hexdigest(),
        "generationConfig": config["generationConfig"], "usageMetadata": result.get("usageMetadata"),
        "finishReason": candidates[0].get("finishReason") if candidates else None})
    return parsed


def extract_article(article, config=None, caller=generate):
    config = config or json.loads((HERE / "model_settings.json").read_text(encoding="utf-8"))
    prompts = json.loads((HERE / "prompts.json").read_text(encoding="utf-8"))
    text = '"""\n' + format_article_content(article) + '\n"""'
    audit = []
    coarse = caller(prompts["coarse_extraction_model"], text, config, audit, "coarse_extraction")
    template = prompts["coarse_refinement_user_template"]
    refinement = template.replace("{{PRE_EXTRACTED_DATA}}", json.dumps(coarse, ensure_ascii=False)).replace("{{ARTICLE}}", text)
    refined_coarse = caller(prompts["coarse_extraction_refinement_model"], refinement, config, audit, "coarse_refinement")
    reactions = extract_reactions(refined_coarse)
    reactions = reactions[:3] if len(reactions) >= 6 else reactions
    details = []
    for i, reaction in enumerate(reactions, 1):
        mode = reaction.get("reactionType")
        prefix = {"thermal catalysis": "thermal", "others": "thermal", "photocatalysis": "photo", "electrocatalysis": "electro"}.get(mode)
        if prefix is None:
            continue
        specified = json.dumps(reaction, ensure_ascii=False)
        detail_user = f"Specified reaction:\n{specified}\nBased on the information in the following paper, you need to extract the information corresponding to the above Specified reaction and complete the code of the reaction:\n{text}"
        extracted = caller(prompts[prefix + "_detail_extraction_model"], detail_user, config, audit, f"detail_extraction_{i}")
        refine_user = f"Specified reaction:\n{json.dumps(extracted, ensure_ascii=False)}\nBased on the information in the following paper, you need to extract the information corresponding to the above specified reaction and complete the missing details of the reaction. In particular, check for missing conversion yield selectivity values.\n{text}"
        refined = caller(prompts[prefix + "_detail_refine_model"], refine_user, config, audit, f"detail_refinement_{i}")
        details.append({"ID": i, "refinedReaction": refined})
    info = {k: v for k, v in article.items() if k in ("title", "doi", "journalName", "pubDate", "citations", "impactEffect")}
    return {"articleInformation": info, "mainReactions": refined_coarse, "reactionDetails": details}, {
        "created_at_UTC": datetime.now(timezone.utc).isoformat(), "api": config["api"],
        "article_input_sha256": hashlib.sha256(json.dumps(article, sort_keys=True).encode("utf-8")).hexdigest(),
        "calls": audit}


def validate(root, directory, output):
    root, directory, output = Path(root), Path(directory), Path(output)
    config = json.loads(Path(root / 'extraction/model_settings.json').read_text(encoding='utf-8'))
    prompts = json.loads(Path(root / 'extraction/prompts.json').read_text(encoding='utf-8'))
    with (root / 'data/samples.csv').open(encoding='utf-8', newline='') as handle:
        samples = list(csv.DictReader(handle))
    if len(samples) != 40 or len({row['doi'].lower() for row in samples}) != 40:
        raise RuntimeError('The sample must contain 40 distinct articles')

    checks = []
    for sample in samples:
        path = directory / Path(sample['article_file']).name
        article = json.loads(path.read_text(encoding='utf-8'))
        if str(article.get('doi', '')).lower() != sample['doi'].lower():
            raise RuntimeError(f'DOI mismatch: {path.name}')
        if not article.get('paragraphs'):
            raise RuntimeError(f'No article paragraphs: {path.name}')
        checksum = hashlib.sha256(path.read_bytes()).hexdigest()
        if checksum != sample['sha256']:
            raise RuntimeError(f'Article checksum mismatch: {path.name}')
        text = format_article_content(article)
        request = build_request(prompts['coarse_extraction_model'], '"""\n' + text + '\n"""', config)
        if request['generationConfig']['responseMimeType'] != 'application/json':
            raise RuntimeError('Extraction requests must use JSON responses')
        checks.append({'sample_id': sample['sample_id'], 'doi': sample['doi'],
                       'file': path.name, 'sha256': checksum, 'paragraphs': len(article['paragraphs'])})

    result = {'articles': len(checks), 'model': config['model'],
              'live_model_called': False, 'inputs': checks}
    output.mkdir(parents=True, exist_ok=True)
    (output / 'article_checks.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--live", action="store_true", help="Explicitly send the supplied article to Gemini")
    args = parser.parse_args()
    if not args.live:
        parser.error("Add --live to run Gemini extraction")
    article = json.loads(args.input.read_text(encoding="utf-8"))
    extracted, audit = extract_article(article)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(extracted, ensure_ascii=False, indent=2), encoding="utf-8")
    args.output.with_suffix(".audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved extraction and audit under {args.output.parent}")
