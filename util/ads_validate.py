#!/usr/bin/env python3
"""Check the hand-written AAS bibliographies against citation keys and ADS.

The bibliographies stay as written: this script never overwrites
paper/references.bib or paper/references-seti.bib. A full run writes
paper/bibcodes.json (citation key to ADS bibcode) and a comparison dump at
paper/ads_export_raw.bib.

    python util/ads_validate.py --keys    # citation keys resolve; no network
    python util/ads_validate.py           # also check each entry against ADS

ADS_TOKEN is read from the environment or from .env at the repo root.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIBS = [ROOT / "paper" / "references.bib", ROOT / "paper" / "references-seti.bib"]
TEX = ROOT / "paper" / "main.tex"
BIBCODES_PATH = ROOT / "paper" / "bibcodes.json"
EXPORT_PATH = ROOT / "paper" / "ads_export_raw.bib"
API = "https://api.adsabs.harvard.edu/v1"

CITE = re.compile(r"\\(?:no)?cite[a-zA-Z]*\*?(?:\s*\[[^\]]*\]){0,2}\s*\{([^}]*)\}")
COMMENT = re.compile(r"(?<!\\)%.*")
ENTRY_START = re.compile(r"@(\w+)\s*\{")

SKIP_KINDS = {"comment", "string", "preamble"}


def load_env(path: Path) -> None:
    """Set variables from a .env file without overriding the environment."""
    if not path.is_file():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ.setdefault(key.strip(), value.strip().strip("\"'"))


def parse_fields(body: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    i = 0
    n = len(body)
    while i < n:
        while i < n and body[i] in " \t\r\n,":
            i += 1
        if i >= n:
            break
        name_match = re.match(r"(\w+)\s*=\s*", body[i:])
        if not name_match:
            i += 1
            continue
        name = name_match.group(1).lower()
        i += name_match.end()
        if i >= n:
            break
        if body[i] == "{":
            depth = 0
            start = i + 1
            while i < n:
                if body[i] == "{":
                    depth += 1
                elif body[i] == "}":
                    depth -= 1
                    if depth == 0:
                        i += 1
                        break
                i += 1
            fields[name] = body[start : i - 1].strip()
        elif body[i] == '"':
            i += 1
            start = i
            while i < n and body[i] != '"':
                i += 1
            fields[name] = body[start:i].strip()
            i += 1
        else:
            start = i
            while i < n and body[i] not in ",\n":
                i += 1
            fields[name] = body[start:i].strip()
    return fields


def parse_bib(path: Path) -> dict[str, tuple[str, dict[str, str]]]:
    text = path.read_text(encoding="utf-8")
    entries: dict[str, tuple[str, dict[str, str]]] = {}
    i = 0
    while True:
        at = text.find("@", i)
        if at < 0:
            break
        match = ENTRY_START.match(text, at)
        if not match:
            i = at + 1
            continue
        kind = match.group(1).lower()
        body_start = match.end()
        if kind in SKIP_KINDS:
            i = body_start
            continue
        comma = text.find(",", body_start)
        if comma < 0:
            break
        key = text[body_start:comma].strip()
        depth = 1
        j = comma + 1
        while j < len(text) and depth:
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
            j += 1
        entries[key] = (kind, parse_fields(text[comma + 1 : j - 1]))
        i = j
    return entries


def load_bibs() -> tuple[dict[str, tuple[str, dict[str, str], str]], list[str]]:
    entries: dict[str, tuple[str, dict[str, str], str]] = {}
    duplicates: list[str] = []
    for path in BIBS:
        if not path.is_file():
            raise SystemExit(f"missing {path.relative_to(ROOT)}")
        for key, (kind, fields) in parse_bib(path).items():
            if key in entries:
                duplicates.append(key)
            entries[key] = (kind, fields, path.name)
    return entries, duplicates


def cite_keys() -> list[str]:
    text = COMMENT.sub("", TEX.read_text(encoding="utf-8"))
    keys: list[str] = []
    seen: set[str] = set()
    for group in CITE.findall(text):
        for key in group.split(","):
            key = key.strip()
            if key and key not in seen:
                seen.add(key)
                keys.append(key)
    return keys


def norm(value: str) -> str:
    value = re.sub(r"\\[a-zA-Z]+", "", value)
    value = re.sub(r"[{}$~\\]", "", value).lower()
    return re.sub(r"[^a-z0-9 ]", "", value)


def ads_request(method: str, path: str, payload: dict | None = None, params: dict | None = None):
    token = os.environ.get("ADS_TOKEN", "")
    if not token:
        raise SystemExit(
            "ADS_TOKEN is not set. Put it in .env at the repo root "
            "(https://ui.adsabs.harvard.edu/user/settings/token)."
        )
    url = API + path
    data = None
    headers = {"Authorization": "Bearer " + token}
    if params:
        url += "?" + urllib.parse.urlencode(params)
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", "replace")[:400]
        raise SystemExit(f"ADS returned {exc.code} for {path}: {detail}") from exc


def search(query: str) -> list[dict]:
    fl = "bibcode,title,author,year,pub,volume,page,doi,identifier,doctype"
    result = ads_request("GET", "/search/query", params={"q": query, "fl": fl, "rows": 5})
    return result["response"]["docs"]


def doi_of(fields: dict[str, str]) -> str:
    doi = fields.get("doi", "")
    return re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi, flags=re.I)


def find(fields: dict[str, str]) -> tuple[list[dict], str]:
    doi = doi_of(fields)
    if doi:
        docs = search(f'doi:"{doi}"')
        if docs:
            return docs, "doi"
    if "eprint" in fields:
        eprint = re.sub(r"(?i)^arxiv:", "", fields["eprint"]).strip()
        docs = search(f"arXiv:{eprint}")
        if docs:
            return docs, "arxiv"
    author = re.sub(r"[{}]", "", fields.get("author", "")).split(" and ")[0]
    author = author.split(",")[0].strip()
    year = fields.get("year", "")
    title = norm(fields.get("title", ""))
    words = " ".join(title.split()[:8])
    if author and year and words:
        docs = search(f'author:"^{author}" year:{year} title:"{words}"')
        if docs:
            return docs, "title"
    if author and year:
        lo = int(year) - 1
        hi = int(year) + 1
        loose = " ".join(title.split()[:4])
        docs = search(f'author:"^{author}" year:{lo}-{hi} abs:"{loose}"')
        return docs, "loose"
    return [], "unsearchable"


def issues_for(fields: dict[str, str], doc: dict) -> list[str]:
    issues = []
    if fields.get("year") and str(doc.get("year")) != fields["year"]:
        issues.append(f"year {fields['year']}->{doc.get('year')}")
    if fields.get("volume") and str(doc.get("volume")) != fields["volume"]:
        issues.append(f"vol {fields['volume']}->{doc.get('volume')}")
    page = (doc.get("page") or [""])[0]
    if fields.get("pages") and norm(fields["pages"].split("-")[0]) != norm(page):
        issues.append(f"page {fields['pages']}->{page}")
    ads_title = doc["title"][0] if doc.get("title") else ""
    if norm(fields.get("title", ""))[:40] != norm(ads_title)[:40]:
        issues.append(f'TITLE: "{ads_title}"')
    first = norm(re.sub(r"[{}]", "", fields.get("author", "")).split(",")[0])
    ads_author = doc["author"][0] if doc.get("author") else ""
    if first and first not in norm(ads_author):
        issues.append(f"AUTHOR {first} -> {ads_author}")
    ads_doi = (doc.get("doi") or [""])[0]
    if ads_doi and doi_of(fields) and ads_doi.lower() != doi_of(fields).lower():
        issues.append(f"doi -> {ads_doi}")
    return issues


def check_keys(entries: dict, duplicates: list[str]) -> list[str]:
    problems = []
    if duplicates:
        problems.append("duplicate keys: " + ", ".join(duplicates))
    missing = [key for key in cite_keys() if key not in entries]
    if missing:
        problems.append("cited but absent from the bibs: " + ", ".join(missing))
    print(f"bib entries: {len(entries)}")
    print(f"citation keys: {len(cite_keys())}")
    for problem in problems:
        print(problem)
    return problems


def validate_ads(entries: dict) -> list[str]:
    report = []
    bibcodes: dict[str, str] = {}
    problems = []
    for key, (_kind, fields, _source) in entries.items():
        docs, how = find(fields)
        time.sleep(0.1)
        if not docs:
            title = fields.get("title", "")
            report.append((key, "NOT FOUND", how, title))
            problems.append(key)
            continue
        docs.sort(key=lambda doc: (doc.get("doctype") == "eprint",))
        doc = docs[0]
        issues = issues_for(fields, doc)
        report.append((key, doc["bibcode"], how, "; ".join(issues) if issues else "OK"))
        bibcodes[key] = doc["bibcode"]
        if issues:
            problems.append(key)
    for key, bibcode, how, message in report:
        print(f"{key:<28} {bibcode:<22} {how:<12} {message}")
    BIBCODES_PATH.write_text(
        json.dumps(dict(sorted(bibcodes.items())), indent=2) + "\n", encoding="utf-8"
    )
    if bibcodes:
        exported = ads_request(
            "POST",
            "/export/bibtex",
            payload={
                "bibcode": list(bibcodes.values()),
                "keyformat": "%1H:%Y",
                "maxauthor": 0,
                "journalformat": 1,
            },
        )
        header = (
            "% Raw ADS BibTeX for comparison only.\n"
            "% Do not replace paper/references.bib or paper/references-seti.bib\n"
            "% with this file: those keep the AAS journal macros.\n\n"
        )
        EXPORT_PATH.write_text(header + exported["export"], encoding="utf-8")
    print(f"exported {len(bibcodes)}; problems {len(problems)}")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keys",
        action="store_true",
        help="check citation keys against the bibs and do not contact ADS",
    )
    args = parser.parse_args()
    load_env(ROOT / ".env")
    entries, duplicates = load_bibs()
    problems = check_keys(entries, duplicates)
    if args.keys:
        return 1 if problems else 0
    problems.extend(validate_ads(entries))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
