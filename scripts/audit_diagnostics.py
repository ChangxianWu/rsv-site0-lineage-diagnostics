"""Corrected RSV diagnostic-component reproduction from supplied windows and metadata.

python -I audit_diagnostics.py --input-dir INPUTS --output-dir NEW_RESULTS
python -I audit_diagnostics.py --self-test-only

Inputs are read-only. This component does not retrieve data or reproduce the entire study.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime, timezone
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import time

import numpy as np
import regex

METADATA = "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv"
TARGET_YEAR = 2026
START_YEAR = 2015
DRIFT_START = 2018
DRIFT_MIN_N = 30
MONTHS = {m.lower(): i for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1)}
MISSING = {"", "na", "n/a", "nan", "none", "null", "not collected", "not provided",
           "missing", "unknown", "not available", "-"}
COMP = str.maketrans("ACGT", "TGCA")
ARMS = {
    "M": {
        "F3": "GCTGTTCAATACAATGTCCTAGA", "B3": "GGTAAATTTGCTGGGCATT",
        "FIP_F1c": "TCTGCTGGCATGGATGATTGG", "FIP_F2": "AGACGATGATCCTGCATCA",
        "BIP_B1c": "CTAGTGAAACAAATATCCACACCCA", "BIP_B2": "GCACTGCACTTCTTGAGTT",
        "LF": "ACATGGGCACCCATATTGTAAG",
    },
    "L": {
        "F3": "TTTGGTGGTGGTGATCCT", "B3": "CAAACTCGGCATTTGGATT",
        "FIP_F1c": "AGCTCAACACAAACACTGAATGTAC", "FIP_F2": "TATATCGAAGCTTTTATAGGAGAAC",
        "BIP_B1c": "CTGGTCACGATCTACAAGATAAGC", "BIP_B2": "TGTGATGACACAAGTCAAGAA",
        "LF": "TATAGCTTCTGTAAGGAAGTCTGGA", "LB": "GGATCTTCCAGATGATAGACTGAAC",
    },
    "N": {
        "Old_N_5331890": "AAGATGCAAATCATAAATTCACAGGA",
        "New_N_5331890": "GCATCATTGCTGTCATTAACTCAATT",
    },
}
COHORTS = {
    "M": ("win_M_full.fasta", "A", True),
    "L": ("win_L_B.fasta", "B", False),
    "N": ("win_n_full.fasta", "B", True),
}
POSITION_KEYS = [f"frac_pos_-{i}" for i in range(1, 6)] + ["frac_pos_internal", "frac_pos_5p"]


def rc(sequence):
    return sequence.translate(COMP)[::-1]


def parse_date_year(raw):
    """Return (year, precision/format, status); never infer a year from a suffix.

    Same-year ranges resolve a year. Cross-year ranges have no unique year and
    are unparseable for this annual analysis; publication fallback is explicit.
    """
    text = str(raw or "").strip()
    if text.lower() in MISSING:
        return None, "missing", "missing"
    if "/" in text:
        parts = text.split("/")
        if len(parts) != 2:
            return None, "unsupported_range", "unparseable"
        left, right = [parse_date_year(s) for s in parts]
        if left[0] is not None and right[0] is not None:
            if left[0] == right[0]:
                return left[0], "same_year_range", "parsed"
            return None, "cross_year_range", "unparseable"
        return None, "invalid_range", "unparseable"
    matches = [
        (r"(\d{4})", "year", "y"),
        (r"(\d{4})-(\d{2})-(\d{2})", "ISO_day", "ymd"),
        (r"(\d{4})-(\d{2})", "ISO_month", "ym"),
        (r"(\d{4})-([A-Za-z]{3})", "year_Mon", "yb"),
        (r"([A-Za-z]{3})-(\d{4})", "Mon_year", "by"),
        (r"(\d{1,2})-([A-Za-z]{3})-(\d{4})", "day_Mon_year", "dby"),
    ]
    for pattern, label, order in matches:
        match = re.fullmatch(pattern, text)
        if not match:
            continue
        values = dict(zip(order, match.groups()))
        try:
            year = int(values["y"])
            month = MONTHS[values["b"].lower()] if "b" in values else int(values.get("m", 1))
            day = int(values.get("d", 1))
            date(year, month, day)
        except (KeyError, ValueError):
            return None, "invalid_" + label, "unparseable"
        return year, label, "parsed"
    return None, "unsupported_format", "unparseable"


def canonical_date(row):
    cy, cf, cs = parse_date_year(row.get("collection_date", ""))
    py, pf, ps = parse_date_year(row.get("first_public", ""))
    if cy is not None:
        year, source, reason = cy, "collection_date", ""
    elif py is not None:
        year, source, reason = py, "first_public_fallback", f"collection_{cs}:{cf}"
    else:
        year, source, reason = None, "unresolved", f"collection_{cs}:{cf};public_{ps}:{pf}"
    return {
        "canonical_year": year, "year_source": source, "fallback_reason": reason,
        "collection_year": cy, "collection_date_format": cf, "collection_date_status": cs,
        "first_public_year": py, "first_public_format": pf,
        "collection_date_raw": row.get("collection_date", ""),
        "first_public_raw": row.get("first_public", ""),
    }


@dataclass(frozen=True)
class Hit:
    orientation: str
    start: int  # zero-based, half-open source coordinates
    end: int
    substring: str
    mismatch_indices: tuple[int, ...]  # zero-based in oriented query
    positions_3prime: tuple[int, ...]  # 1 = terminal, in original oligo orientation

    @property
    def distance(self):
        return len(self.mismatch_indices)


def make_hit(query, substring, orientation, start):
    assert len(query) == len(substring)
    indices = tuple(i for i in range(len(query)) if query[i] != substring[i])
    positions = tuple(sorted(len(query) - i if orientation == "forward" else i + 1
                             for i in indices))
    return Hit(orientation, start, start + len(query), substring, indices, positions)


class Matcher:
    """All globally minimum canonical substitution-only hits at distance <=2.

    Overlapped iteration retains tied sites and both orientations. Hamming
    distances are independently checked after the regex candidate search.
    """
    def __init__(self, oligo):
        if not oligo or set(oligo) - set("ACGT"):
            raise ValueError("Oligos must contain only A/C/G/T")
        self.oligo = oligo
        self.queries = [("forward", oligo), ("reverse", rc(oligo))]
        self.patterns = [(orientation, query, regex.compile(
            f"(?:{query}){{s<=2,i<=0,d<=0}}")) for orientation, query in self.queries]

    def search(self, window):
        exact = []
        for orientation, query in self.queries:
            start = window.find(query)
            while start != -1:
                exact.append(make_hit(query, query, orientation, start))
                start = window.find(query, start + 1)
        if exact:
            return tuple(sorted(exact, key=lambda h: (h.orientation != "forward", h.start)))
        best_distance, hits = 3, []
        for orientation, query, pattern in self.patterns:
            for match in pattern.finditer(window, overlapped=True):
                segment = match.group()
                if set(segment) - set("ACGT"):
                    continue
                hit = make_hit(query, segment, orientation, match.start())
                assert match.fuzzy_counts == (hit.distance, 0, 0)
                assert hit.end == match.end() and hit.distance <= 2
                if hit.distance < best_distance:
                    best_distance, hits = hit.distance, [hit]
                elif hit.distance == best_distance:
                    hits.append(hit)
        return tuple(sorted(hits, key=lambda h: (h.orientation != "forward", h.start)))


def hit_flags(hit, length):
    pos = set(hit.positions_3prime)
    return {
        "p_3prox": any(d <= 5 for d in pos),
        **{f"frac_pos_-{d}": d in pos for d in range(1, 6)},
        "frac_pos_internal": any(5 < d <= length - 5 for d in pos),
        "frac_pos_5p": any(d > length - 5 for d in pos),
    }


def summarize_hits(hits, length):
    if not hits:
        return {"matched": False, "distance": None, "tie_position_disagreement": False,
                "flags_any": {k: False for k in ["p_3prox"] + POSITION_KEYS},
                "flags_all": {k: False for k in ["p_3prox"] + POSITION_KEYS}}
    flags = [hit_flags(h, length) for h in hits]
    return {
        "matched": True, "distance": hits[0].distance,
        "tie_position_disagreement": len({h.positions_3prime for h in hits}) > 1,
        "flags_any": {k: any(f[k] for f in flags) for k in flags[0]},
        "flags_all": {k: all(f[k] for f in flags) for k in flags[0]},
    }


def brute_search(window, oligo):
    hits = []
    for orientation, query in [("forward", oligo), ("reverse", rc(oligo))]:
        for start in range(len(window) - len(query) + 1):
            segment = window[start:start + len(query)]
            if set(segment) - set("ACGT"):
                continue
            hit = make_hit(query, segment, orientation, start)
            if hit.distance <= 2:
                hits.append(hit)
    if not hits:
        return ()
    minimum = min(h.distance for h in hits)
    return tuple(sorted((h for h in hits if h.distance == minimum),
                        key=lambda h: (h.orientation != "forward", h.start)))


def run_self_tests():
    passed = []
    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        passed.append(name)

    for raw, expected in [("2026-01-31", 2026), ("2026", 2026), ("2026-Oct", 2026),
                          ("Oct-2026", 2026), ("2026-10", 2026), ("03-Jan-2020", 2020),
                          ("2020-02-29", 2020), ("2019-02-29", None),
                          ("2015/2016", None), ("2020/2020", 2020),
                          ("Nov-2022/Feb-2023", None), ("", None), ("nonsense2026", None)]:
        check("date:" + raw, parse_date_year(raw)[0] == expected)
    check("collection_precedes_publication", canonical_date({"collection_date": "2019-01-03",
          "first_public": "2026-01-01"})["canonical_year"] == 2019)
    check("cross_year_range_explicit_fallback", canonical_date({"collection_date": "2015/2016",
          "first_public": "2026-01-01"})["year_source"] == "first_public_fallback")
    check("unresolved_without_either_date", canonical_date({})["canonical_year"] is None)
    oligo = "ACGTGCTAGTCACTGATCGGATCA"
    matcher = Matcher(oligo)
    def mutate(seq, indices):
        seq = list(seq)
        for i in indices:
            seq[i] = next(b for b in "ACGT" if b != seq[i])
        return "".join(seq)
    for n in range(4):
        altered = mutate(oligo, [4, 10, 17][:n])
        hits = matcher.search(altered)
        check(f"exactly_{n}_substitutions", bool(hits) == (n <= 2)
              and (not hits or hits[0].distance == n))
        check(f"hamming_reference_{n}", hits == brute_search(altered, oligo))
    for orientation in ["forward", "reverse"]:
        query = oligo if orientation == "forward" else rc(oligo)
        for terminal_position in [1, 2, 3, 4, 5, 6, len(oligo)]:
            index = len(query) - terminal_position if orientation == "forward" else terminal_position - 1
            hits = matcher.search(mutate(query, [index]))
            check(f"{orientation}_position_{terminal_position}", len(hits) == 1
                  and hits[0].orientation == orientation
                  and hits[0].positions_3prime == (terminal_position,))
    two_internal = summarize_hits(matcher.search(mutate(oligo, [8, 10])), len(oligo))
    check("internal_is_record_boolean_not_two_events", two_internal["flags_any"]["frac_pos_internal"] is True)
    two_5p = summarize_hits(matcher.search(mutate(oligo, [1, 3])), len(oligo))
    check("5p_is_record_boolean_not_two_events", two_5p["flags_any"]["frac_pos_5p"] is True)
    for label, seq in [("insertion", oligo[:10] + "A" + oligo[10:]),
                       ("deletion", oligo[:10] + oligo[11:]),
                       ("ambiguous", oligo[:10] + "N" + oligo[11:])]:
        hits = matcher.search(seq)
        check(label + "_not_falsely_exact", not hits or hits[0].distance > 0)
        check(label + "_reference", hits == brute_search(seq, oligo))
    check("ambiguous_candidate_excluded", not matcher.search(oligo[:10] + "N" + oligo[11:]))
    check("canonical_hit_allowed_in_ambiguous_window", matcher.search("NN" + oligo + "NN")[0].distance == 0)
    tied_window = mutate(oligo, [len(oligo)-1]) + "NNNN" + mutate(oligo, [8])
    tied = matcher.search(tied_window)
    flags = summarize_hits(tied, len(oligo))
    check("all_tied_best_sites_retained", len(tied) == 2 and tied == brute_search(tied_window, oligo))
    check("ties_report_position_uncertainty", flags["tie_position_disagreement"]
          and flags["flags_any"]["p_3prox"] and not flags["flags_all"]["p_3prox"])
    rng = np.random.default_rng(20261119)
    for trial in range(80):
        query = "".join(rng.choice(list("ACGT"), size=18))
        embedded = mutate(query, rng.choice(18, size=trial % 4, replace=False))
        if trial % 2:
            embedded = rc(embedded)
        window = "".join(rng.choice(list("ACGTN"), size=23)) + embedded + "".join(
            rng.choice(list("ACGTN"), size=19))
        check(f"random_exhaustive_reference_{trial}", Matcher(query).search(window) == brute_search(window, query))
    return {"passed": len(passed), "failed": 0, "tests": passed}


def read_fasta(path):
    header, pieces = None, []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(pieces).upper()
                header, pieces = line[1:].strip(), []
            elif line.strip():
                if header is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                pieces.append(line.strip())
    if header is not None:
        yield header, "".join(pieces).upper()


def load_metadata(path):
    result = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"accession", "collection_date", "first_public"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Metadata requires columns: {', '.join(sorted(required))}")
        for row in reader:
            accession = row["accession"].split(".")[0]
            if not accession or accession in result:
                raise ValueError(f"Empty or duplicate metadata accession: {accession}")
            result[accession] = canonical_date(row)
    return result


def source_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def write_tsv(path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def add_aggregate(counts, summary, ambiguous, n_acgt_starts, length):
    counts["n"] += 1
    counts["n_window_contains_ambiguity"] += ambiguous
    counts["n_no_canonical_candidate_substring"] += n_acgt_starts == 0
    counts["n_window_too_short"] += length == 0
    if summary["matched"]:
        distance = summary["distance"]
        counts["n_exact"] += distance == 0
        counts["n_le1"] += distance <= 1
        counts["n_le2"] += 1
        counts["n_tied_position_disagreement"] += summary["tie_position_disagreement"]
        for key, value in summary["flags_any"].items():
            counts["count_" + key] += value
        for key, value in summary["flags_all"].items():
            counts["count_all_" + key] += value
    else:
        counts["n_no_canonical_hit_le2"] += 1


def aggregate_row(key, counts):
    basis, gene, oligo, subtype, year = key
    n = counts["n"]
    result = {"year_basis": basis, "gene": gene, "oligo": oligo, "subtype": subtype, "year": year,
              "n": n, "n_collection_date": counts["n_collection_date"],
              "n_first_public_fallback": counts["n_first_public_fallback"],
              "n_exact": counts["n_exact"], "n_le1": counts["n_le1"], "n_le2": counts["n_le2"],
              "p_exact": counts["n_exact"] / n, "p_le1": counts["n_le1"] / n,
              "p_le2": counts["n_le2"] / n}
    for key in ["p_3prox"] + POSITION_KEYS:
        result[key] = counts["count_" + key] / n
        result[key + "_all_best_hits"] = counts["count_all_" + key] / n
    for key in ["n_window_contains_ambiguity", "n_no_canonical_candidate_substring", "n_window_too_short",
                "n_no_canonical_hit_le2", "n_tied_position_disagreement"]:
        result[key] = counts[key]
    assert 0 <= result["p_exact"] <= result["p_le1"] <= result["p_le2"] <= 1
    assert result["p_3prox"] <= result["p_le2"]
    assert counts["n_le2"] + counts["n_no_canonical_hit_le2"] == n
    return result


def drift_rows(yearly_rows):
    groups = defaultdict(list)
    for row in yearly_rows:
        if DRIFT_START <= row["year"] <= TARGET_YEAR:
            groups[(row["year_basis"], row["gene"], row["oligo"], row["subtype"])].append(row)
    output = []
    for key, rows in sorted(groups.items()):
        kept = sorted((r for r in rows if r["n"] >= DRIFT_MIN_N), key=lambda r: r["year"])
        years = [r["year"] for r in kept]
        eligible = len(years) >= 3 and max(years) - min(years) >= 3
        slope, intercept = None, None
        if eligible:
            x = np.array(years, dtype=float)
            y = np.array([r["p_exact"] * 100 for r in kept])
            weights = np.array([r["n"] for r in kept], dtype=float)
            center = np.average(x, weights=weights)
            ycenter = np.average(y, weights=weights)
            slope = float(np.sum(weights * (x-center) * (y-ycenter)) / np.sum(weights * (x-center)**2))
            intercept = float(ycenter - slope * center)
        lookup = {r["year"]: r for r in rows}
        output.append({
            "year_basis": key[0], "gene": key[1], "oligo": key[2], "subtype": key[3],
            "status": "descriptive_WLS" if eligible else "insufficient_years_or_span",
            "slope_pp_per_calendar_year": slope, "intercept_percent": intercept,
            "n_years_used": len(years), "years_used": ";".join(map(str, years)),
            "min_records_per_year": DRIFT_MIN_N,
            "n_used": sum(r["n"] for r in kept),
            "n_first_public_fallback_used": sum(r["n_first_public_fallback"] for r in kept),
            "p_exact_2018": lookup.get(2018, {}).get("p_exact"),
            "p_exact_2026": lookup.get(2026, {}).get("p_exact"),
        })
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, help="Directory containing the four required input files")
    parser.add_argument("--output-dir", type=Path, help="New or empty result directory, outside inputs and release")
    parser.add_argument("--self-test-only", "--selftest", action="store_true", help="Run synthetic checks without input data or output files")
    args = parser.parse_args()
    if not args.self_test_only and (args.input_dir is None or args.output_dir is None):
        parser.error("--input-dir and --output-dir are required unless --self-test-only is used")
    tests = run_self_tests()
    print(f"Synthetic validation: {tests['passed']} passed, 0 failed", flush=True)
    if args.self_test_only:
        return
    input_dir, output_dir = args.input_dir.resolve(), args.output_dir.resolve()
    release_dir = Path(__file__).resolve().parent
    if not input_dir.is_dir():
        parser.error(f"Input directory does not exist: {input_dir}")
    if output_dir == input_dir or input_dir in output_dir.parents or output_dir in input_dir.parents:
        parser.error("Input and output directories must be separate and not nested")
    if output_dir == release_dir or release_dir in output_dir.parents:
        parser.error("Keep generated data outside the code release directory")
    inputs = [input_dir / METADATA] + [input_dir / c[0] for c in COHORTS.values()]
    missing = [str(p) for p in inputs if not p.is_file()]
    if missing:
        parser.error("Missing required input files: " + ", ".join(missing))
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        parser.error("Output directory must be new or empty; existing results are not overwritten")
    started, started_utc = time.monotonic(), datetime.now(timezone.utc).isoformat()
    before = [source_hash(p) for p in inputs]
    metadata = load_metadata(input_dir / METADATA)
    output_dir.mkdir(parents=True, exist_ok=True)
    date_counts = Counter()
    def count_dates(scope, gene, info):
        date_counts[(scope, gene, info["year_source"], info["canonical_year"],
                     info["collection_date_status"], info["collection_date_format"],
                     info["first_public_format"], info["fallback_reason"])] += 1
    for info in metadata.values():
        count_dates("all_metadata", "", info)
    totals = defaultdict(Counter)
    cohort_provenance, cohort_dates = [], []
    audit_fields = ["gene", "oligo", "subtype", "accession", "canonical_year", "year_source",
                    "fallback_reason", "collection_year", "window_length", "window_non_acgt_count",
                    "canonical_substring_starts", "min_hamming_le2", "status", "best_hit_count",
                    "tie_position_disagreement", "selected_orientation", "selected_start_0based",
                    "selected_end_exclusive", "selected_substring", "selected_3prime_positions_1based",
                    "best_hits_json"] + ["p_3prox"] + POSITION_KEYS
    with gzip.open(output_dir / "per_record_audit.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene, (filename, subtype, require_tag) in COHORTS.items():
            print(f"Reading {gene}/{subtype} cohort from {filename}", flush=True)
            records, seen, stats = [], set(), Counter()
            for header, sequence in read_fasta(input_dir / filename):
                stats["input_records"] += 1
                parts = header.split()
                tags = dict(p.split("=", 1) for p in parts[1:] if "=" in p)
                tag = tags.get("subtype", parts[1] if len(parts) > 1 and "=" not in parts[1] else "")
                if require_tag and not tag:
                    raise ValueError(f"Missing subtype: {header}")
                if tag and tag != subtype:
                    stats["excluded_other_subtype"] += 1
                    continue
                accession = parts[0].split(".")[0]
                if accession in seen:
                    raise ValueError(f"Duplicate cohort accession: {gene}/{accession}")
                seen.add(accession)
                stats["selected_subtype_records"] += 1
                info = metadata.get(accession)
                if info is None:
                    info = canonical_date({})
                    info["fallback_reason"] = "metadata_accession_missing"
                    stats["metadata_accession_missing"] += 1
                count_dates("selected_subtype_windows", gene, info)
                cohort_dates.append({"gene": gene, "subtype": subtype, "accession": parts[0], **info})
                year = info["canonical_year"]
                if year is None or not START_YEAR <= year <= TARGET_YEAR:
                    stats["excluded_unresolved_or_outside_2015_2026"] += 1
                    continue
                stats["analysis_records_2015_2026"] += 1
                count_dates("analysis_2015_2026", gene, info)
                if year == TARGET_YEAR:
                    count_dates("analysis_2026", gene, info)
                records.append((parts[0], sequence, info))
            cohort_provenance.append({"gene": gene, "subtype": subtype, "file": filename,
                "subtype_selection": "explicit FASTA subtype tag" if require_tag else
                "supplied subtype-B-only L file; checked against any available tags", **dict(stats)})
            for oligo, sequence in ARMS[gene].items():
                matcher, cache = Matcher(sequence), {}
                arm_started = time.monotonic()
                for accession, window, info in records:
                    cached = cache.get(window)
                    if cached is None:
                        hits = matcher.search(window)
                        summary = summarize_hits(hits, len(sequence))
                        ambiguous = sum(base not in "ACGT" for base in window)
                        canonical_starts = sum(max(0, len(m.group()) - len(sequence) + 1)
                                               for m in re.finditer("[ACGT]+", window))
                        cached = (hits, summary, ambiguous, canonical_starts)
                        cache[window] = cached
                    hits, summary, ambiguous, canonical_starts = cached
                    for basis in ["canonical_collection_first"] + (
                            ["collection_date_only"] if info["year_source"] == "collection_date" else []):
                        key = (basis, gene, oligo, subtype, info["canonical_year"])
                        add_aggregate(totals[key], summary, bool(ambiguous), canonical_starts,
                                      0 if len(window) < len(sequence) else len(window))
                        totals[key]["n_" + ("collection_date" if info["year_source"] == "collection_date"
                                           else "first_public_fallback")] += 1
                    selected = hits[0] if hits else None
                    status = "canonical_hit_le2" if hits else (
                        "no_full_length_canonical_substring" if canonical_starts == 0 else "no_canonical_hit_le2")
                    writer.writerow({
                        "gene": gene, "oligo": oligo, "subtype": subtype, "accession": accession,
                        "canonical_year": info["canonical_year"], "year_source": info["year_source"],
                        "fallback_reason": info["fallback_reason"], "collection_year": info["collection_year"],
                        "window_length": len(window), "window_non_acgt_count": ambiguous,
                        "canonical_substring_starts": canonical_starts,
                        "min_hamming_le2": summary["distance"], "status": status,
                        "best_hit_count": len(hits), "tie_position_disagreement": int(summary["tie_position_disagreement"]),
                        "selected_orientation": selected.orientation if selected else "",
                        "selected_start_0based": selected.start if selected else "",
                        "selected_end_exclusive": selected.end if selected else "",
                        "selected_substring": selected.substring if selected else "",
                        "selected_3prime_positions_1based": ";".join(map(str, selected.positions_3prime)) if selected else "",
                        "best_hits_json": json.dumps([{
                            "orientation": h.orientation, "start_0based": h.start,
                            "end_exclusive": h.end, "substring": h.substring,
                            "mismatch_positions_3prime_1based": h.positions_3prime,
                        } for h in hits], separators=(",", ":")),
                        **{k: int(v) for k, v in summary["flags_any"].items()},
                    })
                print(f"  {oligo}: {len(records):,} records, {len(cache):,} unique windows; "
                      f"{time.monotonic()-arm_started:.1f}s", flush=True)
    yearly = [aggregate_row(key, counts) for key, counts in sorted(totals.items())]
    denominator_sets = defaultdict(set)
    for row in yearly:
        denominator_sets[(row["year_basis"], row["gene"], row["year"])].add(row["n"])
    assert all(len(values) == 1 for values in denominator_sets.values())
    final = [r for r in yearly if r["year"] == TARGET_YEAR and r["year_basis"] == "canonical_collection_first"]
    sensitivity = [r for r in yearly if r["year"] == TARGET_YEAR and r["year_basis"] == "collection_date_only"]
    dates = [{"scope": k[0], "gene": k[1], "year_source": k[2], "canonical_year": k[3],
              "collection_date_status": k[4], "collection_date_format": k[5],
              "first_public_format": k[6], "fallback_reason": k[7], "n_records": n}
             for k, n in sorted(date_counts.items(), key=lambda pair: tuple(str(v) for v in pair[0]))]
    write_tsv(output_dir / "per_arm_2026.tsv", final)
    write_tsv(output_dir / "per_arm_2026_collection_only.tsv", sensitivity)
    write_tsv(output_dir / "per_arm_by_year.tsv", yearly)
    write_tsv(output_dir / "exact_drift_2018_2026.tsv", drift_rows(yearly))
    write_tsv(output_dir / "date_source_counts.tsv", dates)
    write_tsv(output_dir / "cohort_record_dates.tsv", cohort_dates)
    write_json(output_dir / "synthetic_tests.json", tests)
    write_json(output_dir / "methods_and_limitations.json", {
        "cohorts": cohort_provenance, "target_year": TARGET_YEAR, "oligos": ARMS,
        "year_definition": "Collection year first; publication year only for missing/unparseable collection dates. Cross-year ranges do not resolve a collection year. Fallback is NOT collection year; collection-only sensitivity supplied.",
        "denominator": "Fixed distinct-accession gene/subtype/year cohort; ambiguous, short and unmatched windows remain in n. Missing gene windows are outside the cohort.",
        "search": "Both orientations, all overlapping canonical substrings, substitutions only {s<=2,i<=0,d<=0}; all globally minimum hits <=2 retained. Non-ACGT substrings excluded.",
        "positions": "One-based 3-prime positions: forward L-i; reverse i+1. All -1..-5 included. Internal excludes both terminal five-base regions; internal/5p are per-record Boolean marginals.",
        "ties": "Primary position fractions: ANY tied best hit; *_all_best_hits: ALL tied best hits. Representative hit is deterministic; all hits remain in audit.",
        "position_denominator": "Same fixed n; positions observed only in canonical best hits <=2. Unmatched contributes no observed position, not proof of no mismatch.",
        "drift": "Descriptive WLS exact percent versus calendar year, 2018-2026, weights n; annual n>=30; >=3 years spanning >=3 years. Collection-only analysis separate; no causal inference or uncertainty model.",
        "limitations": [
            "Component reproduction only: no retrieval, window extraction, subtype reassignment, exposure models, phylogeny or figures.",
            "Existing orientation-normalized windows and supplied subtype assignments are input assumptions; untagged L file must already be subtype B only.",
            "No canonical hit <=2 can reflect ambiguity, missing sites, indels or divergence; it does not establish >=3 biological substitutions.",
            "Whole-window best-site matching does not establish the unique biological binding site.",
            "Sequence compatibility is not assay performance; composite-primer arm roles qualify positional interpretation.",
            "Publication fallback and changing sequencing composition can affect annual trends.",
        ],
    })
    after = [source_hash(p) for p in inputs]
    if before != after:
        raise RuntimeError("Input files changed during analysis")
    write_json(output_dir / "provenance.json", {
        "started_utc": started_utc, "elapsed_seconds": time.monotonic()-started,
        "script": source_hash(Path(__file__)), "python": sys.version,
        "numpy": np.__version__, "regex": regex.__version__,
        "inputs_unchanged": True, "inputs": before,
        "outputs": [source_hash(p) for p in sorted(output_dir.iterdir()) if p.name != "provenance.json"],
    })
    print("\n2026 collection-first results (publication fallback explicitly counted):", flush=True)
    for row in final:
        print(f"{row['gene']} {row['oligo']:16} n={row['n']:5d} fallback={row['n_first_public_fallback']:4d} "
              f"exact={100*row['p_exact']:7.3f}% <=1={100*row['p_le1']:7.3f}% "
              f"<=2={100*row['p_le2']:7.3f}% 3prox={100*row['p_3prox']:6.3f}%", flush=True)
    print(f"Output: {output_dir}\nInputs unchanged; elapsed {time.monotonic()-started:.1f}s", flush=True)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
