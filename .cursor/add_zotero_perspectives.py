#!/usr/bin/env python3
"""Add missing perspectives papers to Zotero (one leaf collection each)."""

from __future__ import annotations

import json
import re
import sqlite3
import time
import urllib.parse
import urllib.request
from pathlib import Path

LIBRARY_ID = Path.home().joinpath(".zotc/library").read_text().strip()
API_KEY = Path.home().joinpath(".zotc/api-key").read_text().strip()
API = f"https://api.zotero.org/users/{LIBRARY_ID}/items"
HEADERS_JSON = {
    "Zotero-API-Key": API_KEY,
    "Content-Type": "application/json",
}

STOP = {
    "a",
    "an",
    "the",
    "of",
    "and",
    "for",
    "in",
    "on",
    "to",
    "with",
    "from",
    "by",
    "at",
    "or",
    "as",
    "via",
    "using",
    "its",
    "into",
    "over",
}

# doi, collection key, notes
TARGETS: list[tuple[str, str]] = [
    # CTS / joint inference
    ("10.1038/s43018-022-00356-3", "LGR4XKTI"),  # BayesPrism
    ("10.1038/s41467-021-26328-2", "LGR4XKTI"),  # BLADE
    ("10.1186/s13619-026-00299-5", "LGR4XKTI"),  # BayesPrism-DWLS
    # alternative distributions
    ("10.1111/biom.13614", "QWC7SK8H"),  # IsoDeconvMM
    ("10.1101/2023.12.07.570524", "QWC7SK8H"),  # DeconV
    ("10.3389/fevo.2021.588292", "QWC7SK8H"),  # PLN versatile
    ("10.1007/s11222-025-10729-0", "QWC7SK8H"),  # ZIPLN
    # covariates / WLS / GLS
    ("10.1093/bib/bbac430", "SU4ZK7S7"),  # MuSiC2
    ("10.1186/s13059-025-03776-3", "SU4ZK7S7"),  # Unico
    ("10.3389/fgene.2022.825896", "SU4ZK7S7"),  # DecOT
    # semi-reference
    ("10.1214/20-AOAS1376", "MZLIESEB"),  # BayICE
    ("10.1158/0008-5472.CAN-25-4102", "MZLIESEB"),  # CDState
    # ensemble
    ("10.1093/bioinformatics/btac825", "NVYXUP5L"),  # EnDecon
    ("10.1038/s41467-024-50618-0", "NVYXUP5L"),  # White DREAM
    # hierarchical / archetypes / potency
    ("10.1002/advs.202514073", "S52P99M7"),  # HIDF
    ("10.1093/bib/bbaf556", "S52P99M7"),  # scDETECT
    ("10.1093/bioinformatics/btae709", "S52P99M7"),  # tissueResolver
    ("10.1038/s41467-020-18416-6", "5X49IBGS"),  # ACTION
    ("10.1371/journal.pcbi.1010025", "5X49IBGS"),  # scAAnet
    ("10.1073/pnas.2530194123", "5X49IBGS"),  # Pareto archetypes 2025
    ("10.1038/s42003-024-07109-1", "5PZ7PFGU"),  # organoid potency
    # dynamics
    ("10.1038/s41564-025-01983-z", "D7KQ3TKA"),  # ChronoStrain
    ("10.64898/2026.01.26.701686", "D7KQ3TKA"),  # scTREND
    # spatial
    ("10.3390/ijms27125259", "I6XYZB4P"),  # HistoMap
    ("10.64898/2026.01.06.697922", "I6XYZB4P"),  # HEDeST
    ("10.64898/2026.02.10.705204", "I6XYZB4P"),  # SpaDecoder
    ("10.1101/2025.08.12.669903", "7VT8VP9K"),  # ST cell-type mismatch
    # multi-omics
    ("10.1038/s41592-026-03007-y", "VTA7RHPX"),  # DECODE
    ("10.1038/s41587-026-03218-w", "VTA7RHPX"),  # ATAC+RNA
    ("10.64898/2026.06.01.729268", "VTA7RHPX"),  # proteomics
    ("10.3390/biotech15010007", "VTA7RHPX"),  # multi-scale TIME
]


def existing_dois() -> set[str]:
    db = Path.home() / "Zotero" / "zotero.sqlite"
    uri = f"file:{db}?mode=ro&immutable=1"
    con = sqlite3.connect(uri, uri=True)
    cur = con.cursor()
    rows = cur.execute(
        """
        SELECT v.value
        FROM itemData d
        JOIN itemDataValues v ON v.valueID = d.valueID
        WHERE d.fieldID = 59
        """
    ).fetchall()
    out = set()
    for (val,) in rows:
        if not val:
            continue
        d = val.strip().lower().replace("https://doi.org/", "").replace("http://doi.org/", "")
        out.add(d)
    return out


def http_json(url: str, timeout: int = 30):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "DeCovarT-zotero-import/1.0 (mailto:bastien_chassagnol@laposte.net)"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as exc:  # noqa: BLE001
        print(f"  fetch fail {url}: {exc}")
        return None


def parse_authors(raw) -> list[dict]:
    creators = []
    if not raw:
        return creators
    for a in raw:
        if isinstance(a, str):
            parts = a.replace(",", " ").split()
            if not parts:
                continue
            if len(parts) == 1:
                creators.append({"creatorType": "author", "lastName": parts[0]})
            else:
                creators.append(
                    {
                        "creatorType": "author",
                        "firstName": " ".join(parts[:-1]),
                        "lastName": parts[-1],
                    }
                )
            continue
        family = a.get("family") or a.get("lastName") or ""
        given = a.get("given") or a.get("firstName") or ""
        if a.get("name") and not family:
            family = a["name"]
        if family:
            creators.append(
                {
                    "creatorType": "author",
                    "firstName": given,
                    "lastName": family,
                }
            )
    return creators


def citekey(creators, title: str, year: str) -> str:
    last = "anon"
    if creators:
        last = re.sub(r"[^A-Za-z-]", "", creators[0].get("lastName", "anon"))
        last = last.lower()
    words = re.findall(r"[A-Za-z0-9]+", title or "")
    kept = []
    for w in words:
        if w.lower() in STOP and kept:
            continue
        kept.append(w[:1].upper() + w[1:].lower())
        if len(kept) >= 4:
            break
    short = "".join(kept) if kept else "Untitled"
    year = re.sub(r"[^0-9]", "", year)[:4] or "nd"
    return f"{last}{short}{year}"


def from_crossref(doi: str) -> dict | None:
    data = http_json(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}")
    if not data or "message" not in data:
        return None
    m = data["message"]
    title = " ".join(m.get("title") or [""]).strip()
    authors = []
    for a in m.get("author") or []:
        authors.append(
            {
                "creatorType": "author",
                "firstName": a.get("given", ""),
                "lastName": a.get("family") or a.get("name", ""),
            }
        )
    year = ""
    for key in ("published-print", "published-online", "issued", "created"):
        parts = (m.get(key) or {}).get("date-parts") or []
        if parts and parts[0]:
            year = str(parts[0][0])
            break
    container = " ".join(m.get("container-title") or [""]).strip()
    item_type = "journalArticle"
    if "posted-content" in (m.get("type") or "") or not container:
        if any(x in (container + title).lower() for x in ("biorxiv", "arxiv")):
            item_type = "preprint"
    item = {
        "itemType": item_type,
        "title": title,
        "creators": authors,
        "date": year,
        "DOI": doi,
        "url": f"https://doi.org/{doi}",
    }
    if container:
        if item_type == "journalArticle":
            item["publicationTitle"] = container
        else:
            item["repository"] = container
    if m.get("volume"):
        item["volume"] = str(m["volume"])
    if m.get("issue"):
        item["issue"] = str(m["issue"])
    if m.get("page"):
        item["pages"] = str(m["page"])
    issn = m.get("ISSN") or []
    if issn:
        item["ISSN"] = issn[0]
    return item


def from_openalex(doi: str) -> dict | None:
    data = http_json(f"https://api.openalex.org/works/doi:{doi}")
    if not data:
        return None
    title = data.get("title") or data.get("display_name") or ""
    authors = []
    for aa in data.get("authorships") or []:
        a = aa.get("author") or {}
        name = a.get("display_name") or ""
        parts = name.split()
        if not parts:
            continue
        authors.append(
            {
                "creatorType": "author",
                "firstName": " ".join(parts[:-1]),
                "lastName": parts[-1],
            }
        )
    year = str(data.get("publication_year") or "")
    loc = data.get("primary_location") or {}
    source = (loc.get("source") or {}).get("display_name") or ""
    item_type = "journalArticle" if source and "rxiv" not in source.lower() else "preprint"
    item = {
        "itemType": item_type,
        "title": title,
        "creators": authors,
        "date": year,
        "DOI": doi,
        "url": f"https://doi.org/{doi}",
    }
    if source:
        if item_type == "journalArticle":
            item["publicationTitle"] = source
        else:
            item["repository"] = source
    bibl = data.get("biblio") or {}
    if bibl.get("volume"):
        item["volume"] = str(bibl["volume"])
    if bibl.get("issue"):
        item["issue"] = str(bibl["issue"])
    if bibl.get("first_page") and bibl.get("last_page"):
        item["pages"] = f"{bibl['first_page']}-{bibl['last_page']}"
    return item if title else None


def from_biorxiv(doi: str) -> dict | None:
    data = http_json(f"https://api.biorxiv.org/details/biorxiv/{doi}")
    if not data:
        data = http_json(f"https://api.biorxiv.org/details/medrxiv/{doi}")
    coll = (data or {}).get("collection") or []
    if not coll:
        return None
    rec = coll[-1]
    authors = []
    for name in (rec.get("authors") or "").split(";"):
        name = name.strip()
        if not name:
            continue
        parts = name.replace(",", " ").split()
        if len(parts) == 1:
            authors.append({"creatorType": "author", "lastName": parts[0]})
        else:
            authors.append(
                {
                    "creatorType": "author",
                    "firstName": " ".join(parts[:-1]),
                    "lastName": parts[-1],
                }
            )
    year = (rec.get("date") or "")[:4]
    return {
        "itemType": "preprint",
        "title": rec.get("title") or "",
        "creators": authors,
        "date": year,
        "DOI": doi,
        "url": f"https://doi.org/{doi}",
        "repository": rec.get("server") or "bioRxiv",
        "abstractNote": rec.get("abstract") or "",
    }


def enrich_item(item: dict) -> dict:
    year = str(item.get("date") or "")[:4]
    ck = citekey(item.get("creators") or [], item.get("title") or "", year)
    extra = item.get("extra") or ""
    lines = [ln for ln in extra.splitlines() if ln.strip()]
    if not any(ln.lower().startswith("citation key:") for ln in lines):
        lines.append(f"Citation Key: {ck}")
    item["extra"] = "\n".join(lines)
    item["citationKey"] = ck
    return item


def post_items(items: list[dict]) -> list[dict]:
    body = json.dumps(items).encode("utf-8")
    req = urllib.request.Request(
        API,
        data=body,
        headers=HEADERS_JSON,
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    return payload


def main() -> None:
    have = existing_dois()
    created = []
    skipped = []
    failed = []
    batch = []
    meta = []

    for doi, coll in TARGETS:
        dlow = doi.lower()
        if dlow in have or any(dlow in h or h in dlow for h in have):
            skipped.append(doi)
            print(f"SKIP already in library: {doi}")
            continue
        print(f"RESOLVE {doi}")
        item = from_crossref(doi) or from_openalex(doi) or from_biorxiv(doi)
        if not item or not item.get("title"):
            failed.append(doi)
            print(f"  NO METADATA {doi}")
            continue
        item = enrich_item(item)
        item["collections"] = [coll]
        print(f"  -> {item.get('citationKey')} | {item.get('title')[:80]}")
        batch.append(item)
        meta.append((doi, item.get("citationKey"), item.get("title")))
        time.sleep(0.2)

    print(f"\nTo create: {len(batch)}; skip {len(skipped)}; fail {len(failed)}")
    # Zotero max 50 per request
    for i in range(0, len(batch), 50):
        chunk = batch[i : i + 50]
        try:
            resp = post_items(chunk)
        except Exception as exc:  # noqa: BLE001
            print("POST failed:", exc)
            # retry one-by-one
            for it in chunk:
                try:
                    resp_one = post_items([it])
                    print("  one-ok", it.get("citationKey"), resp_one)
                    created.append(it)
                except Exception as exc2:  # noqa: BLE001
                    print("  one-fail", it.get("DOI"), exc2)
                    failed.append(it.get("DOI"))
            continue
        print(json.dumps({k: resp.get(k) for k in ("success", "failed", "unchanged") if k in resp}, indent=2))
        if resp.get("failed"):
            print("failed detail", resp["failed"])
        created.extend(chunk)

    out = Path(__file__).with_name("zotero-perspectives-added.json")
    out.write_text(
        json.dumps(
            {
                "created": meta,
                "skipped": skipped,
                "failed": failed,
            },
            indent=2,
        )
        + "\n"
    )
    print("wrote", out)


if __name__ == "__main__":
    main()
