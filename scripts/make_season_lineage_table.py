"""Season-level share of subtype B records carrying the B.D.E.1 marker set.

Recomputes, from frozen inputs, the series plotted in Figure 2B and published as
Supplementary Table S11:

  input 1  lineage_marker_assignment_B_2026-10-28.tsv (frozen; n = 10,174 = marker-callable
           subtype B records with a collection year in 2015-2026)
  input 2  ena_rsv_meta_matrix_sp11250_2026-10-27.tsv (frozen; collection_date month/year)
  rules    season_of is copied verbatim from seasonize_and_h1base.py (pre-registered P7 v0.9):
           northern hemisphere Oct-Dec -> "y/y+1", Jan-May -> "y-1/y", Jun-Sep interseason;
           southern hemisphere Jan-Sep -> "y", Oct-Dec -> "y+1"; tropical countries follow the
           northern rules.

Records with year-only collection dates (no month) and interseason records are excluded from the
series; their counts are printed. Output: lineage_share_by_season_<date>.tsv with season,
convention, n, k_bde1, pct.
"""
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
OUT = WORK / "lineage_share_by_season_2026-11-16.tsv"

MNTH = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06",
        "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}
SOUTH = {"AUS", "ARG", "CHL", "NZL", "ZAF", "URY", "PRY", "BOL"}
TROPICAL = {"KEN", "MOZ", "GMB", "PHL", "VNM", "IND", "SGP", "MYS", "IDN",
            "COL", "PER", "BRA", "ECU", "PAN", "CRI", "NIC", "HND", "GTM"}


def season_of(iso3, coll):
    """Season assignment, verbatim from seasonize_and_h1base.py (no first-public fallback)."""
    month = year = None
    if len(coll) >= 8 and "-" in coll:
        parts = coll.split("-")
        if len(parts) == 3:
            if parts[0].isdigit() and len(parts[0]) == 4:
                year, month = parts[0], MNTH.get(parts[1], parts[1] if parts[1].isdigit() else None)
            else:
                month = MNTH.get(parts[1])
                year = parts[2] if parts[2].isdigit() and len(parts[2]) == 4 else (
                    "20" + parts[2] if len(parts[2]) == 2 else None)
    if month is None or year is None:
        return None
    y, m = int(year), int(month)
    if iso3 in SOUTH and iso3 not in TROPICAL:
        return f"{y}" if m <= 9 else f"{y + 1}"
    if 10 <= m <= 12:
        return f"{y}/{y + 1}"
    if 1 <= m <= 5:
        return f"{y - 1}/{y}"
    return "interseason"


lin = pd.read_csv(WORK / f"lineage_marker_assignment_B_{DAY}.tsv", sep="\t", dtype=str)
assert len(lin) == 10174, f"marker-callable set changed: {len(lin)}"
lin["acc"] = lin.accession.str.split(".").str[0]
lin["iso3"] = lin.iso3.fillna("")
lin["bde1"] = lin.BDE1 == "True"

mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["acc"] = mx.accession.astype(str).str.split(".").str[0]
mx["coll"] = mx.collection_date.fillna("")
mx = mx.sort_values("coll", ascending=False).drop_duplicates("acc", keep="first")
coll_of = mx.set_index("acc")["coll"].to_dict()

agg = defaultdict(lambda: [0, 0])   # (season, convention) -> [k, n]
skipped = {"no_month": 0, "interseason": 0, "no_meta": 0}
for r in lin.itertuples():
    coll = coll_of.get(r.acc, "")
    if coll == "":
        skipped["no_meta"] += 1
        continue
    s = season_of(r.iso3, coll)
    if s is None:
        skipped["no_month"] += 1
        continue
    if s == "interseason":
        skipped["interseason"] += 1
        continue
    conv = "southern" if "/" not in s else "northern"
    a = agg[(s, conv)]
    a[1] += 1
    a[0] += int(r.bde1)

rows = []
for (s, conv), (k, n) in sorted(agg.items()):
    rows.append({"season": s, "convention": conv, "n": n, "k_bde1": k,
                 "pct": round(k / n * 100, 1) if n else ""})
out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False)
print(f"excluded: {skipped}  | seasons: {len(out)}")
north = out[(out.convention == "northern") &
            (out.season.astype(str).str.split("/").str[0].astype(int).between(2015, 2025))]
print(north.to_string(index=False))
