"""Assemble the frozen analysis base (v1.0) from the seasonised base (v0.9).

The v0.9 base carries the per-record season, country and key-mutation flags. This step adds the
four fields the exposure models consume:

  season_start  first year of the collection season (integer)
  y             frozen outcome: subtype B K68Q / K68N / K201S; subtype A K65R / K209E
  lin_b         lineage-composition covariate: I206M carriage in subtype B records
  study         ENA study accession, used as an alternative clustering unit

Output: h1_analysis_base_v1.0_2026-10-28.tsv.gz

This script restates the step that produced the frozen base; the downstream estimates are
reproduced when it is run on the released v0.9 base.
"""
from pathlib import Path

import pandas as pd

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
OUTCOME = {"B": {"K68Q", "K68N", "K201S"}, "A": {"K65R", "K209E"}}

base = pd.read_csv(WORK / f"h1_analysis_base_v0.9_{DAY}.tsv.gz", sep="\t", dtype=str)
mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["accession_base"] = mx["accession"].astype(str).str.split(".").str[0]
study = mx.set_index("accession_base")["study_accession"].fillna("NA").to_dict()

base["accession_base"] = base["accession"].astype(str).str.split(".").str[0]
base["study"] = base["accession_base"].map(study).fillna("NA")
base["season_start"] = base["season"].astype(str).str.split("/").str[0]

y, lin = [], []
for _, r in base.iterrows():
    st = r["subtype"]
    muts = {c[2:] for c in base.columns if c.startswith("m_") and str(r[c]) == "1"}
    y.append(int(bool(OUTCOME[st] & muts)))
    lin.append(int("I206M" in muts) if st == "B" else 0)
base["y"] = y
base["lin_b"] = lin

base.to_csv(WORK / f"h1_analysis_base_v1.0_{DAY}.tsv.gz", sep="\t", index=False,
            compression="gzip")
print("rows:", len(base), "| events:", int(base.y.sum()),
      "| studies:", base.study.nunique(), "| countries:", base.iso3.nunique())
print(base[["iso3", "season_start", "y", "lin_b", "study"]].head().to_string(index=False))
