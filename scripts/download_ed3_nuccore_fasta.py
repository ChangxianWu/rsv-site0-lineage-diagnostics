"""Bulk retrieval of RSV nucleotide records (NCBI Nucleotide, taxid 11250).

Downloads the nuccore FASTA for the human RSV taxon into raw/ (read-only layer), resuming
between batches, and writes a manifest with the file sha256.

Safety: requests are restricted to https + an allow-listed host, redirects are refused, and the
output path is boundary-checked. Sequence orientation is handled downstream (the ENA metadata
matrix supplies mol_type as a prior and the window extractor scores both strands), so no
multi-gigabyte GenBank download is performed here.
"""
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlparse

import requests

WORK = Path(__file__).resolve().parent
RAW = WORK.parent / "raw"
RAW.mkdir(exist_ok=True)
DAY = "2026-10-28"
ES = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ALLOWED_HOSTS = {"eutils.ncbi.nlm.nih.gov"}


def guarded_get(url, **kw):
    p = urlparse(url)
    if p.scheme != "https" or (p.hostname or "").lower() not in ALLOWED_HOSTS:
        raise ValueError(f"blocked non-allowlisted URL: {url}")
    kw["allow_redirects"] = False
    return requests.get(url, **kw)


def out_path(name):
    p = (WORK / name).resolve()
    if p.parent != WORK:
        raise ValueError(f"path escape blocked: {name}")
    return p


RAW.mkdir(exist_ok=True)
OUT_FASTA = RAW / f"h_rsv_nuccore_txid11250_{DAY}.fasta"
IDS_FILE = out_path(f"nuccore_uids_txid11250_{DAY}.txt")
MANIFEST = out_path(f"download_manifest_{DAY}.json")

# 1) UID list (resume: skip esearch when the file already exists)
if IDS_FILE.exists() and IDS_FILE.stat().st_size > 100:
    uids = IDS_FILE.read_text().split()
    print(f"resume: {len(uids)} uids from file", flush=True)
else:
    r = guarded_get(ES + "esearch.fcgi", params={
        "db": "nuccore", "term": "txid11250[Organism]",
        "retmode": "json", "retmax": 90000}, timeout=120)
    uids = r.json()["esearchresult"]["idlist"]
    IDS_FILE.write_text("\n".join(uids))
    print(f"esearch: {len(uids)} uids saved", flush=True)

# 2) Batched efetch with resume (batch count derived from records already written)
BATCH = 200
start = 0
tmp = OUT_FASTA.with_suffix(".part")
if tmp.exists():
    n_have = tmp.read_text().count("\n>")
    start = (n_have // BATCH) * BATCH if n_have else 0
    print(f"resume fasta from batch {start // BATCH} ({n_have} records already)", flush=True)

t0 = time.time()
n_ok = 0
with tmp.open("a", encoding="utf-8") as fh:
    for i in range(start, len(uids), BATCH):
        chunk = uids[i:i + BATCH]
        for attempt in range(5):
            try:
                r = guarded_get(ES + "efetch.fcgi", params={
                    "db": "nuccore", "id": ",".join(chunk),
                    "rettype": "fasta", "retmode": "text"}, timeout=180)
                if r.status_code == 200 and "LOCUS" not in r.text[:200]:
                    txt = r.text if r.text.endswith("\n") else r.text + "\n"
                    fh.write(txt)
                    n_ok += 1
                    break
                raise RuntimeError(f"batch {i}: http {r.status_code}")
            except (requests.RequestException, RuntimeError) as e:
                print(f"batch {i} attempt {attempt + 1}: {e}", flush=True)
                time.sleep(2 * (attempt + 1))
        else:
            print(f"batch {i}: FAILED after retries, stopping for resume", flush=True)
            break
        if (i // BATCH) % 20 == 0:
            print(f"batch {i // BATCH}: {i + len(chunk)}/{len(uids)} uids, "
                  f"{time.time() - t0:.0f}s", flush=True)
        time.sleep(0.45)  # NCBI rate-limit discipline

# 3) Completion check + manifest + sha256
n_records = tmp.read_text().count("\n>") if tmp.exists() else 0
done = n_records >= len(uids) * 0.999
if done:
    tmp.replace(OUT_FASTA)
    h = hashlib.sha256(OUT_FASTA.read_bytes()).hexdigest()
    manifest = {
        "file": OUT_FASTA.name, "n_uids": len(uids), "n_records": n_records,
        "sha256": h, "source": "eutils efetch nuccore rettype=fasta",
        "term": "txid11250[Organism]", "date": DAY, "status": "complete",
    }
else:
    manifest = {
        "file": tmp.name, "n_uids": len(uids), "n_records": n_records,
        "source": "eutils efetch nuccore rettype=fasta",
        "term": "txid11250[Organism]", "date": DAY, "status": "incomplete-resume",
    }
MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=1))
print(json.dumps(manifest, ensure_ascii=False), flush=True)
print("DONE", flush=True)
