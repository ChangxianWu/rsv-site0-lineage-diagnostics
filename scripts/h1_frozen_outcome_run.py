"""H1 run on the frozen outcome set.

Outcome set = the site catalog (v1.0) rows flagged as interim H1 outcomes:
  subtype B: K68Q / K68N / K201S; subtype A: K65R / K209E.

Model: record-level binomial GEE with country-clustered robust standard errors, a linear season
trend and subtype/tropical covariates; sensitivity analyses drop one country at a time and drop
the Kenyan 2016/17 cluster.

Output: h1_frozen_outcome_run_2026-10-28.tsv
"""
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
SOUTH = {"AUS", "ARG", "CHL", "NZL", "ZAF", "URY", "PRY", "BOL"}
TROPICAL = {"KEN", "MOZ", "GMB", "PHL", "VNM", "IND", "COL", "PER", "BRA"}
EXPOSED_FROM = {"ESP": 2023, "FRA": 2023, "LUX": 2023, "USA": 2023,
                "DEU": 2024, "NLD": 2024, "ITA": 2024, "AUT": 2024, "IRL": 2024}
PARTIAL = {"BEL", "PRT"}
BMARK = "I206M"
MNTH = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06",
        "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}

mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["accession_base"] = mx["accession"].astype(str).str.split(".").str[0]
mx["iso3"] = mx["country"].fillna("").str.split(":").str[0].str.strip().map(
    {"Spain": "ESP", "France": "FRA", "United States": "USA", "USA": "USA",
     "Luxembourg": "LUX", "Austria": "AUT", "Ireland": "IRL", "Germany": "DEU",
     "Netherlands": "NLD", "Belgium": "BEL", "Portugal": "PRT", "Italy": "ITA",
     "United Kingdom": "GBR", "Norway": "NOR", "Denmark": "DNK", "Sweden": "SWE",
     "Finland": "FIN", "Switzerland": "CHE", "Australia": "AUS", "Canada": "CAN",
     "Japan": "JPN", "China": "CHN", "Kenya": "KEN", "South Africa": "ZAF",
     "Brazil": "BRA", "Argentina": "ARG", "Mexico": "MEX", "India": "IND",
     "Gambia": "GMB", "Mozambique": "MOZ", "Philippines": "PHL", "Vietnam": "VNM",
     "Chile": "CHL", "Peru": "PER", "Colombia": "COL", "Israel": "ISR",
     "Turkey": "TUR", "Poland": "POL", "Czechia": "CZE", "Greece": "GRC"}).fillna("")
mx["coll"] = mx["collection_date"].fillna("")
mx["fpy"] = mx["first_public"].fillna("").str[:4]
meta = mx.set_index("accession_base")[["iso3", "coll", "fpy"]].to_dict("index")

def season_of(iso3, coll, fpy):
    parts = coll.split("-")
    month = year = None
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
    hemi = "S" if iso3 in SOUTH else "N"
    if hemi == "N":
        if 10 <= m <= 12:
            return f"{y}/{y + 1}"
        if 1 <= m <= 5:
            return f"{y - 1}/{y}"
        return "interseason"
    return f"{y}" if m <= 9 else f"{y + 1}"

per = pd.read_csv(WORK / f"traj_per_record_site0_2026-10-28.tsv.gz", sep="\t", dtype=str)
recs = []
for _, r in per.iterrows():
    st, acc = r["subtype"], r["accession"]
    base_acc = acc.split(".")[0]
    m = meta.get(base_acc, {})
    iso3 = m.get("iso3", "") or ""
    season = season_of(iso3, m.get("coll", ""), m.get("fpy", ""))
    if not season or season == "interseason":
        continue
    ss = int(season.split("/")[0])
    if ss < 2010:
        continue
    muts = set((r.get("mutations") or "").split(";")) - {""}
    if st == "B":
        y = int(bool({"K68Q", "K68N", "K201S"} & muts))
    else:
        y = int(bool({"K65R", "K209E"} & muts))
    lin = int(BMARK in muts) if st == "B" else 0
    recs.append({"country": iso3, "subtype": st, "season_start": ss, "y": y, "lin_b": lin,
                 "tropical": int(iso3 in TROPICAL)})
base = pd.DataFrame(recs)
base["exposed"] = base["country"].map(EXPOSED_FROM).notna().astype(int)
base["exp_start"] = base["country"].map(EXPOSED_FROM)
base["exp_post"] = ((base["exposed"] == 1) &
                    (base["season_start"] >= base["exp_start"])).astype(int)
base["partial"] = base["country"].isin(PARTIAL).astype(int)
print("records:", len(base), "events:", int(base.y.sum()), flush=True)

def run_gee(d, tag):
    if d["y"].sum() == 0:
        return {"tag": tag, "note": "no events"}
    m = smf.gee("y ~ exp_post + lin_b + C(subtype) + tropical + season_start",
                groups="country", data=d, family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Independence()).fit(maxiter=500)
    return {"tag": tag, "n": len(d), "events": int(d.y.sum()),
            "converged": bool(m.converged),
            "OR": round(float(np.exp(m.params["exp_post"])), 3),
            "CI_lo": round(float(np.exp(m.conf_int().loc["exp_post", 0])), 3),
            "CI_hi": round(float(np.exp(m.conf_int().loc["exp_post", 1])), 3),
            "p": round(float(m.pvalues["exp_post"]), 4)}

rows = [run_gee(base[base["partial"] == 0], "main_frozen_outcome")]
for c in sorted(EXPOSED_FROM):
    rows.append(run_gee(base[(base["partial"] == 0) & (base["country"] != c)], f"loo_{c}"))
rows.append(run_gee(base[(base["partial"] == 0) & (base["country"] != "KEN")], "no_KEN"))
out = pd.DataFrame([{k: v for k, v in r.items()} for r in rows])
out.to_csv(WORK / f"h1_frozen_outcome_run_{DAY}.tsv", sep="\t", index=False)
print(out.to_string(index=False), flush=True)
print("saved (exposure table v0.9; still a pre-freeze dry run)", flush=True)
