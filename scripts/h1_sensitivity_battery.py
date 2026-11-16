"""H1 sensitivity battery (the full specification set reported in Table 3).

Model: logit(y) ~ exposed x post + subtype + tropical + season_start, with the lineage covariate
(lin_b) added in the primary specification. Output h1_battery_2026-10-28.tsv reports the primary
model, the unadjusted-for-lineage model, difference-in-differences with main effects, study-
clustered standard errors, study fixed effects, leave-one-batch-out (LOBO), country and study
cluster bootstrap intervals, exposure-label permutation, partial-as-exposed and the crude 2x2.

Random draws use a fixed seed, so the file is regenerable byte-for-byte in the table values.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
rng = np.random.default_rng(20261028)

TROP = {"KEN", "MOZ", "GMB", "PHL", "VNM", "IND", "COL", "PER", "BRA"}
EXPOSED_FROM = {"ESP": 2023, "FRA": 2023, "LUX": 2023, "USA": 2023,
                "DEU": 2024, "NLD": 2024, "ITA": 2024, "AUT": 2024, "IRL": 2024}
PARTIAL = {"BEL", "PRT"}

b = pd.read_csv(WORK / f"h1_analysis_base_v1.0_{DAY}.tsv.gz", sep="\t")
b["y"] = pd.to_numeric(b["y"], errors="coerce").fillna(0).astype(int)
b["season_start"] = pd.to_numeric(b["season_start"], errors="coerce")
b = b[b.season_start.notna()].copy()
b["country"] = b["iso3"].fillna("")
b["sexposed"] = b["country"].map(EXPOSED_FROM).notna().astype(int)
b["exp_start"] = b["country"].map(EXPOSED_FROM)
b["exp_post"] = ((b.sexposed == 1) & (b.season_start >= b.exp_start)).astype(int)
b["post"] = (b.season_start >= 2023).astype(int)
b["tropical"] = b["country"].isin(TROP).astype(int)
b["study"] = b["study"].fillna("NA")
d = b[~b.country.isin(PARTIAL)].copy()

rows = []


def gee(formula, data, groups, term="exp_post"):
    m = smf.gee(formula, groups=groups, data=data, family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Independence()).fit(maxiter=500)
    lo, hi = np.exp(m.conf_int().loc[term])
    return (round(float(np.exp(m.params[term])), 3), round(float(lo), 3),
            round(float(hi), 3), round(float(m.pvalues[term]), 4))


f_main = "y ~ exp_post + lin_b + C(subtype) + tropical + season_start"
f_nolin = "y ~ exp_post + C(subtype) + tropical + season_start"
rows.append({"spec": "main_country_cluster", **dict(zip(
    ["OR", "lo", "hi", "p"], gee(f_main, d, "country"))), "n": len(d), "events": int(d.y.sum())})

f_did = "y ~ sexposed * post + C(subtype) + tropical + season_start"
rows.append({"spec": "main_without_lineage_covariate", **dict(zip(
    ["OR", "lo", "hi", "p"], gee(f_nolin, d, "country"))), "n": len(d), "events": int(d.y.sum())})

rows.append({"spec": "did_with_main_effects", **dict(zip(
    ["OR", "lo", "hi", "p"], gee(f_did, d, "country", term="sexposed:post"))),
    "n": len(d), "events": int(d.y.sum())})

rows.append({"spec": "study_cluster", **dict(zip(
    ["OR", "lo", "hi", "p"], gee(f_main, d, "study"))), "n": len(d), "events": int(d.y.sum())})

# Study fixed effects (GLM with country-clustered errors)
m = smf.glm(f_main + " + C(study)", data=d, family=sm.families.Binomial()).fit(
    cov_type="cluster", cov_kwds={"groups": d["country"]})
lo, hi = np.exp(m.conf_int().loc["exp_post"])
rows.append({"spec": "study_fixed_effects", "OR": round(float(np.exp(m.params["exp_post"])), 3),
             "lo": round(float(lo), 3), "hi": round(float(hi), 3),
             "p": round(float(m.pvalues["exp_post"]), 4), "n": len(d), "events": int(d.y.sum())})

# Partially exposed countries treated as exposed
d2 = b.copy()
d2["exp_post"] = ((d2.sexposed == 1) | d2.country.isin(PARTIAL)).astype(int) * (
    d2.season_start >= d2.exp_start.fillna(2024)).astype(int)
rows.append({"spec": "partial_as_exposed", **dict(zip(
    ["OR", "lo", "hi", "p"], gee(f_main, d2, "country"))), "n": len(d2),
    "events": int(d2.y.sum())})

# LOBO over studies that contribute events
ev_studies = d[d.y == 1].study.value_counts()
lobo = []
for s in ev_studies.index:
    sub = d[d.study != s]
    r = gee(f_main, sub, "country")
    lobo.append((s, *r))
rows.append({"spec": "LOBO_range",
             "OR": f"{min(x[1] for x in lobo):.2f}-{max(x[1] for x in lobo):.2f}",
             "lo": "", "hi": "", "p": "", "n": len(d), "events": int(d.y.sum())})

# Cluster bootstrap (300 replicates), country and study level
def boot(unit, B=300):
    units = d[unit].unique()
    est = []
    for _ in range(B):
        pick = rng.choice(units, size=len(units), replace=True)
        parts = [d[d[unit] == u] for u in pick]
        dd = pd.concat(parts, ignore_index=True)
        if dd.y.sum() == 0:
            continue
        try:
            est.append(gee(f_main, dd, "country")[0])
        except Exception:
            continue
    est = np.array(est)
    return (len(est), round(float(np.percentile(est, 2.5)), 2),
            round(float(np.percentile(est, 97.5)), 2),
            round(float((est > 2).mean()), 3))

nb, blo, bhi, p2 = boot("country")
rows.append({"spec": "country_bootstrap",
             "lo": blo, "hi": bhi, "p": f"P(OR>2)={p2}", "n": f"B={nb}", "events": int(d.y.sum())})
nb, blo, bhi, p2 = boot("study")
rows.append({"spec": "study_bootstrap", "OR": "", "lo": blo, "hi": bhi,
             "p": f"P(OR>2)={p2}", "n": f"B={nb}", "events": int(d.y.sum())})

# Permutation of exposure labels (300 replicates), permuted at country level
obs = gee(f_main, d, "country")[0]
countries = d.country.unique()
perm_or = []
for _ in range(300):
    lab = dict(zip(countries, rng.permutation(
        [EXPOSED_FROM.get(c, 0) for c in countries])))
    dd = d.copy()
    dd["ps"] = dd.country.map(lab)
    dd["exp_post"] = ((dd.ps > 0) & (dd.season_start >= dd.ps)).astype(int)
    try:
        perm_or.append(gee(f_main, dd, "country")[0])
    except Exception:
        continue
perm_or = np.array(perm_or)
rows.append({"spec": "permutation_labels", "OR": obs, "lo": "", "hi": "",
             "p": f"perm_p={(perm_or >= obs).mean():.3f}; frac>=2={(perm_or >= 2).mean():.3f}",
             "n": len(perm_or), "events": int(d.y.sum())})

# Crude 2x2
cr = b[~b.country.isin(PARTIAL)].groupby("exp_post").y.agg(["sum", "size"])
c0, n0 = cr.loc[0, "sum"], cr.loc[0, "size"]
c1, n1 = cr.loc[1, "sum"], cr.loc[1, "size"]
rows.append({"spec": "crude_2x2", "OR": round((c1 / n1) / (c0 / n0), 3),
             "lo": "", "hi": "", "p": "",
             "n": f"unexp {c0}/{n0}; exp-post {c1}/{n1}", "events": int(c0 + c1)})

out = pd.DataFrame(rows)
out.to_csv(WORK / f"h1_battery_{DAY}.tsv", sep="\t", index=False)
print(out.to_string(index=False))
print("\nLOBO detail:")
for s, OR, lo, hi, p in lobo:
    print(f"  drop {s}: OR {OR} ({lo}-{hi})")
