"""Season assignment and the H1 analysis base table.

Inputs: traj_per_record_site0_2026-10-28.tsv.gz (sparse per-record substitutions) and the ENA
metadata matrix (collection_date, country).

Season rules (pre-registered P7 v0.9):
  Northern hemisphere: Oct-Dec -> "y/y+1"; Jan-May -> "y-1/y"; Jun-Sep -> interseason
  (excluded from seasonal aggregation).
  Southern hemisphere: Jan-Sep -> "y"; Oct-Dec -> "y+1" (the following season); no interseason.
  Tropical or unclassifiable: northern rules, with a tropical/unknown flag retained.

Outputs: traj_key_sites_by_season_2026-10-28.tsv and h1_analysis_base_v0.9_2026-10-28.tsv.gz
(accession, subtype, iso3, hemisphere, year, season, group, key-mutation flags).
"""
from collections import defaultdict, Counter
from pathlib import Path

import pandas as pd

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
ANCHORS = {"ESP", "FRA", "USA", "LUX", "AUT", "IRL", "DEU", "NLD", "BEL", "PRT"}
KEYMUTS = {"A": ["K65R", "E66K", "N88S", "A74T", "K68E", "K68Q", "K68N"],
           "B": ["K68Q", "K68N", "I206M", "Q209R", "S211N", "S190N", "K191R"]}
SOUTH = {"AUS", "ARG", "CHL", "NZL", "ZAF", "URY", "PRY", "BOL"}
TROPICAL = {"KEN", "MOZ", "GMB", "PHL", "VNM", "IND", "SGP", "MYS", "IDN", "COL", "PER", "BRA", "ECU", "PAN", "CRI", "NIC", "HND", "GTM"}

CN = {"Spain": "ESP", "France": "FRA", "United States": "USA", "USA": "USA",
      "Luxembourg": "LUX", "Austria": "AUT", "Ireland": "IRL", "Germany": "DEU",
      "Netherlands": "NLD", "Belgium": "BEL", "Portugal": "PRT", "Italy": "ITA",
      "United Kingdom": "GBR", "Norway": "NOR", "Denmark": "DNK", "Sweden": "SWE",
      "Finland": "FIN", "Switzerland": "CHE", "Australia": "AUS", "Canada": "CAN",
      "Japan": "JPN", "China": "CHN", "Kenya": "KEN", "South Africa": "ZAF",
      "Brazil": "BRA", "Argentina": "ARG", "Mexico": "MEX", "India": "IND",
      "Gambia": "GMB", "Mozambique": "MOZ", "Philippines": "PHL", "Vietnam": "VNM",
      "Chile": "CHL", "Peru": "PER", "Colombia": "COL", "Israel": "ISR",
      "Turkey": "TUR", "Poland": "POL", "Czechia": "CZE", "Greece": "GRC"}

mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["accession_base"] = mx["accession"].astype(str).str.split(".").str[0]
mx["iso3"] = mx["country"].fillna("").str.split(":").str[0].str.strip().map(CN).fillna("")
mx["coll"] = mx["collection_date"].fillna("")
mx["fp_year"] = mx["first_public"].fillna("").str[:4]
meta = mx.set_index("accession_base")[["iso3", "coll", "fp_year"]].to_dict("index")

def season_of(iso3, coll, fp_year):
    ym = coll.split("-")  # GenBank "21-Sep-2021", ISO "2021-09-21", or "2021"
    month = None
    year = None
    if len(coll) >= 8 and "-" in coll:
        parts = coll.split("-")
        if len(parts) == 3:
            if parts[0].isdigit() and len(parts[0]) == 4:
                year, month = parts[0], {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
                    "May": "05", "Jun": "06", "Jul": "07", "Aug": "08", "Sep": "09",
                    "Oct": "10", "Nov": "11", "Dec": "12"}.get(parts[1], parts[1] if parts[1].isdigit() else None)
            else:
                month = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05",
                    "Jun": "06", "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10",
                    "Nov": "11", "Dec": "12"}.get(parts[1])
                year = parts[2] if parts[2].isdigit() and len(parts[2]) == 4 else (
                    ("20" + parts[2]) if len(parts[2]) == 2 else None)
    if month is None or year is None:
        return None, fp_year if fp_year else None
    y = int(year)
    m = int(month)
    hemi = "S" if iso3 in SOUTH else ("N" if iso3 else "N")
    if iso3 in TROPICAL:
        hemi = "N"  # tropical countries follow northern rules; a tropical flag is carried separately
    if hemi == "N":
        if 10 <= m <= 12:
            return f"{y}/{y + 1}", y
        if 1 <= m <= 5:
            return f"{y - 1}/{y}", y
        return "interseason", y
    else:
        if 1 <= m <= 9:
            return f"{y}", y
        return f"{y + 1}", y

per = pd.read_csv(WORK / f"traj_per_record_site0_2026-10-28.tsv.gz", sep="\t", dtype=str)
rows = []
agg = defaultdict(Counter)
den = defaultdict(Counter)
for _, r in per.iterrows():
    st, acc = r["subtype"], r["accession"]
    base = acc.split(".")[0]
    m = meta.get(base, {})
    iso3 = m.get("iso3", r.get("iso3", "")) or ""
    season, cal_year = season_of(iso3, m.get("coll", ""), m.get("fp_year", ""))
    if not season or season == "interseason":
        continue
    group = "anchor" if iso3 in ANCHORS else "other"
    den[(st, season, group)]["den"] += 1
    muts = set((r.get("mutations") or "").split(";")) - {""}
    flags = {k: (1 if k in muts else 0) for k in KEYMUTS[st]}
    rows.append({"accession": acc, "subtype": st, "iso3": iso3, "year": cal_year,
                 "season": season, "group": group, **{f"m_{k}": v for k, v in flags.items()}})
    for k in KEYMUTS[st]:
        if flags[k]:
            agg[(st, k, season, group)]["n"] += 1

base = pd.DataFrame(rows)
base.to_csv(WORK / f"h1_analysis_base_v0.9_{DAY}.tsv.gz", sep="\t", index=False, compression="gzip")
out = []
for (st, mut, season, group), c in agg.items():
    d = den[(st, season, group)]["den"]
    out.append({"subtype": st, "mutation": mut, "season": season, "group": group,
                "n": c["n"], "den": d, "freq": round(c["n"] / d, 4)})
df = pd.DataFrame(out).sort_values(["subtype", "mutation", "group", "season"])
df.to_csv(WORK / f"traj_key_sites_by_season_{DAY}.tsv", sep="\t", index=False, encoding="utf-8-sig")
print("base rows:", len(base), "| season-aggregated rows:", len(df), flush=True)
for st in ("A", "B"):
    for mut in KEYMUTS[st][:4]:
        sub = df[(df.subtype == st) & (df.mutation == mut)]
        if len(sub):
            piv = sub.pivot_table(index="group", columns="season", values="freq")
            print(f"--- {st} {mut} ---")
            print(piv.to_string())
