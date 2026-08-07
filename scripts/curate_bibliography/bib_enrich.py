#!/usr/bin/env python3
"""Enrich incomplete BibTeX entries via Crossref (DOI-first, then title search).

Non-interactive CLI intended for Cursor agents and pre-commit workflows.

Examples:
  scripts/curate_bibliography/bib_enrich/.venv/bin/python scripts/curate_bibliography/bib_enrich.py inst/REFERENCES.bib \\
    --write --backup --min-score 92 --report .cursor/bib-enrichment.json

  # dry-run (report only):
  scripts/curate_bibliography/bib_enrich/.venv/bin/python scripts/curate_bibliography/bib_enrich.py inst/REFERENCES.bib \\
    --report .cursor/bib-enrichment.json
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import bibtexparser
import httpx
from rapidfuzz import fuzz
from rich.console import Console

console = Console(stderr=True)

REQUIRED_BY_TYPE: dict[str, set[str]] = {
    "article": {"title", "author", "year", "journal"},
    "book": {"title", "author", "year", "publisher"},
    "incollection": {"title", "author", "year", "booktitle", "publisher"},
    "inproceedings": {"title", "author", "year", "booktitle"},
    "phdthesis": {"title", "author", "year", "school"},
    "mastersthesis": {"title", "author", "year", "school"},
    "techreport": {"title", "author", "year", "institution"},
    "misc": {"title"},
    "online": {"title"},
}

# Fields we may fill from Crossref (never overwrite existing non-empty values).
ENRICHABLE_FIELDS = {
    "title",
    "author",
    "year",
    "journal",
    "booktitle",
    "publisher",
    "doi",
    "volume",
    "number",
    "pages",
    "isbn",
    "issn",
    "editor",
    "edition",
    "school",
    "institution",
    "url",
}


def norm_doi(raw: str | None) -> str | None:
    if not raw:
        return None
    d = raw.strip().strip("{}").strip('"').rstrip(",").strip()
    d = re.sub(r"^https?://(dx\.)?doi\.org/", "", d, flags=re.I)
    return d.lower() or None


def norm_title(raw: str | None) -> str:
    if not raw:
        return ""
    t = re.sub(r"[{}\\$]", "", raw)
    t = re.sub(r"\s+", " ", t).strip().lower()
    return t


def first_author_surname(author: str | None) -> str:
    if not author:
        return ""
    first = author.split(" and ")[0].strip()
    if "," in first:
        return first.split(",", 1)[0].strip().lower()
    parts = first.split()
    return parts[-1].lower() if parts else ""


def missing_required(entry: dict[str, str]) -> list[str]:
    etype = entry.get("ENTRYTYPE", "misc").lower()
    required = REQUIRED_BY_TYPE.get(etype, {"title"})
    # Books may use editor instead of author.
    if etype == "book" and not entry.get("author") and entry.get("editor"):
        required = (required - {"author"}) | {"editor"}
    missing = [f for f in sorted(required) if not (entry.get(f) or "").strip()]
    return missing


def authors_from_crossref(message: dict[str, Any]) -> str | None:
    authors = message.get("author") or []
    parts: list[str] = []
    for a in authors:
        family = (a.get("family") or "").strip()
        given = (a.get("given") or "").strip()
        if family and given:
            parts.append(f"{family}, {given}")
        elif family:
            parts.append(family)
        elif given:
            parts.append(given)
    return " and ".join(parts) if parts else None


def editors_from_crossref(message: dict[str, Any]) -> str | None:
    editors = message.get("editor") or []
    parts: list[str] = []
    for a in editors:
        family = (a.get("family") or "").strip()
        given = (a.get("given") or "").strip()
        if family and given:
            parts.append(f"{family}, {given}")
        elif family:
            parts.append(family)
    return " and ".join(parts) if parts else None


def year_from_crossref(message: dict[str, Any]) -> str | None:
    for key in ("published-print", "published-online", "issued", "created"):
        parts = (message.get(key) or {}).get("date-parts") or []
        if parts and parts[0]:
            return str(parts[0][0])
    return None


def title_from_crossref(message: dict[str, Any]) -> str | None:
    titles = message.get("title") or []
    return titles[0] if titles else None


def container_from_crossref(message: dict[str, Any]) -> str | None:
    containers = message.get("container-title") or []
    return containers[0] if containers else None


def isbn_from_crossref(message: dict[str, Any]) -> str | None:
    isbns = message.get("ISBN") or []
    return isbns[0] if isbns else None


def pages_from_crossref(message: dict[str, Any]) -> str | None:
    return message.get("page")


def map_crossref_to_fields(message: dict[str, Any], entry_type: str) -> dict[str, str]:
    etype = entry_type.lower()
    out: dict[str, str] = {}
    title = title_from_crossref(message)
    if title:
        out["title"] = title
    authors = authors_from_crossref(message)
    if authors:
        out["author"] = authors
    editors = editors_from_crossref(message)
    if editors:
        out["editor"] = editors
    year = year_from_crossref(message)
    if year:
        out["year"] = year
    container = container_from_crossref(message)
    publisher = message.get("publisher")
    doi = message.get("DOI")
    if doi:
        out["doi"] = str(doi).lower()
    if message.get("volume"):
        out["volume"] = str(message["volume"])
    if message.get("issue"):
        out["number"] = str(message["issue"])
    pages = pages_from_crossref(message)
    if pages:
        out["pages"] = str(pages).replace("-", "--")
    isbn = isbn_from_crossref(message)
    if isbn:
        out["isbn"] = str(isbn)
    if message.get("edition-number"):
        out["edition"] = str(message["edition-number"])
    if publisher:
        out["publisher"] = str(publisher)
        # Tech reports / theses sometimes store org as publisher.
        if etype in {"techreport"}:
            out.setdefault("institution", str(publisher))
        if etype in {"phdthesis", "mastersthesis"}:
            out.setdefault("school", str(publisher))

    if etype == "article" and container:
        out["journal"] = container
    elif etype in {"incollection", "inproceedings", "book"} and container:
        out["booktitle"] = container
    elif etype == "book" and not container and title:
        # monograph: container often absent
        pass
    elif container and etype not in {"article"}:
        out.setdefault("booktitle", container)

    return {k: v for k, v in out.items() if k in ENRICHABLE_FIELDS and v}


class CrossrefClient:
    def __init__(self, mailto: str = "bastien_chassagnol@laposte.net", pause_s: float = 0.2):
        self.client = httpx.Client(
            base_url="https://api.crossref.org",
            headers={"User-Agent": f"DeCovarT-bib-enrich/1.0 (mailto:{mailto})"},
            timeout=30.0,
        )
        self.pause_s = pause_s

    def close(self) -> None:
        self.client.close()

    def _get(self, path: str, params: dict[str, Any] | None = None) -> dict[str, Any] | None:
        time.sleep(self.pause_s)
        try:
            r = self.client.get(path, params=params)
            r.raise_for_status()
            return r.json()
        except Exception as exc:  # noqa: BLE001
            console.print(f"[yellow]Crossref request failed[/]: {path} ({exc})")
            return None

    def by_doi(self, doi: str) -> dict[str, Any] | None:
        payload = self._get(f"/works/{doi}")
        if not payload:
            return None
        return payload.get("message")

    def search(self, query: str, rows: int = 5) -> list[dict[str, Any]]:
        payload = self._get(
            "/works",
            params={"query.bibliographic": query, "rows": rows},
        )
        if not payload:
            return []
        return list((payload.get("message") or {}).get("items") or [])


def score_candidate(entry: dict[str, str], message: dict[str, Any]) -> float:
    entry_title = norm_title(entry.get("title"))
    cand_title = norm_title(title_from_crossref(message) or "")
    if not entry_title or not cand_title:
        return 0.0
    title_score = float(fuzz.token_sort_ratio(entry_title, cand_title))
    # Soft bonuses for author / year agreement.
    bonus = 0.0
    entry_year = (entry.get("year") or "").strip()
    cand_year = year_from_crossref(message) or ""
    if entry_year and cand_year and entry_year == cand_year:
        bonus += 3.0
    entry_author = first_author_surname(entry.get("author") or entry.get("editor"))
    cand_authors = authors_from_crossref(message) or ""
    if entry_author and entry_author in cand_authors.lower():
        bonus += 3.0
    return min(100.0, title_score + bonus)


def choose_candidate(
    entry: dict[str, str],
    client: CrossrefClient,
    min_score: float,
) -> tuple[dict[str, Any] | None, str, float, str]:
    """Return (message, method, score, status)."""
    doi = norm_doi(entry.get("doi"))
    if doi:
        msg = client.by_doi(doi)
        if msg:
            score = score_candidate(entry, msg)
            # DOI hit is trusted even if title fuzzy score is imperfect.
            return msg, "doi", max(score, 99.0), "matched"

    title = (entry.get("title") or "").strip()
    author = first_author_surname(entry.get("author") or entry.get("editor"))
    year = (entry.get("year") or "").strip()
    query = " ".join(x for x in [title, author, year] if x)
    if not query:
        return None, "none", 0.0, "unresolved"

    items = client.search(query, rows=5)
    if not items:
        return None, "search", 0.0, "unresolved"

    scored = [(score_candidate(entry, m), m) for m in items]
    scored.sort(key=lambda x: x[0], reverse=True)
    best_score, best = scored[0]
    second = scored[1][0] if len(scored) > 1 else 0.0
    if best_score < min_score:
        return None, "search", best_score, "unresolved"
    if second >= min_score and (best_score - second) < 3.0:
        return None, "search", best_score, "ambiguous"
    return best, "search", best_score, "matched"


def merge_missing(
    entry: dict[str, str],
    proposed: dict[str, str],
) -> dict[str, str]:
    added: dict[str, str] = {}
    for field, value in proposed.items():
        if field not in ENRICHABLE_FIELDS:
            continue
        current = (entry.get(field) or "").strip()
        if current:
            continue
        entry[field] = value
        added[field] = value
    # Promote howpublished leftovers when booktitle/publisher now present.
    how = (entry.get("howpublished") or "").strip()
    if how and entry.get("booktitle") and entry.get("publisher"):
        # Keep howpublished only if it still adds information; otherwise drop.
        if norm_title(how) in {
            norm_title(entry.get("booktitle")),
            norm_title(f"{entry.get('booktitle')} ({entry.get('edition')})"),
        }:
            entry.pop("howpublished", None)
    return added


def enrich_bib(
    path: Path,
    *,
    write: bool,
    backup: bool,
    min_score: float,
    report_path: Path | None,
    only_incomplete: bool = True,
) -> int:
    raw = path.read_text(encoding="utf-8")
    db = bibtexparser.loads(raw)
    client = CrossrefClient()
    report: dict[str, Any] = {
        "source": str(path),
        "min_score": min_score,
        "entries": [],
        "summary": {
            "scanned": 0,
            "incomplete": 0,
            "enriched": 0,
            "unresolved": 0,
            "ambiguous": 0,
            "skipped_complete": 0,
        },
    }

    try:
        for entry in db.entries:
            report["summary"]["scanned"] += 1
            key = entry.get("ID", "?")
            missing = missing_required(entry)
            if only_incomplete and not missing:
                report["summary"]["skipped_complete"] += 1
                continue
            report["summary"]["incomplete"] += 1
            message, method, score, status = choose_candidate(entry, client, min_score)
            row: dict[str, Any] = {
                "key": key,
                "type": entry.get("ENTRYTYPE"),
                "missing_before": missing,
                "method": method,
                "score": score,
                "status": status,
                "fields_added": {},
            }
            if status == "matched" and message is not None:
                proposed = map_crossref_to_fields(message, entry.get("ENTRYTYPE", "misc"))
                added = merge_missing(entry, proposed)
                row["fields_added"] = added
                still = missing_required(entry)
                row["missing_after"] = still
                if added:
                    report["summary"]["enriched"] += 1
                if still:
                    row["status"] = "partial"
                    report["summary"]["unresolved"] += 1
                else:
                    row["status"] = "enriched"
            elif status == "ambiguous":
                report["summary"]["ambiguous"] += 1
            else:
                report["summary"]["unresolved"] += 1
            report["entries"].append(row)
            console.print(
                f"{key}: {row['status']} "
                f"(missing={missing}, score={score:.1f}, method={method})"
            )
    finally:
        client.close()

    if report_path is not None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        console.print(f"Wrote report: {report_path}")

    if write:
        if backup:
            bak = path.with_suffix(path.suffix + ".bak")
            shutil.copy2(path, bak)
            console.print(f"Backup: {bak}")
        writer = bibtexparser.bwriter.BibTexWriter()
        writer.indent = "  "
        writer.order_entries_by = ("ID",)
        path.write_text(writer.write(db), encoding="utf-8")
        console.print(f"Wrote enriched bibliography: {path}")

    unresolved = report["summary"]["unresolved"] + report["summary"]["ambiguous"]
    return 1 if unresolved else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bib", type=Path, help="Path to .bib file")
    parser.add_argument("--write", action="store_true", help="Write enriched .bib in place")
    parser.add_argument("--backup", action="store_true", help="Write .bib.bak before --write")
    parser.add_argument("--min-score", type=float, default=92.0)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument(
        "--all",
        action="store_true",
        help="Attempt enrichment for all entries, not only incomplete ones",
    )
    args = parser.parse_args(argv)
    if not args.bib.is_file():
        console.print(f"[red]Missing bibliography[/]: {args.bib}")
        return 2
    return enrich_bib(
        args.bib,
        write=args.write,
        backup=args.backup,
        min_score=args.min_score,
        report_path=args.report,
        only_incomplete=not args.all,
    )


if __name__ == "__main__":
    sys.exit(main())
