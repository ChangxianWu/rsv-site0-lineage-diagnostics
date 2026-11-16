"""Pilot audit of the N-gene probe pair (PMC5331890, 2017).

Anchor: both probes match the subtype B reference (NC_001781) exactly (old at 0-based 164, new
at 925) and differ from the subtype A reference by 3-5 bases, which makes them subtype specific.

Audit: for subtype B records, the exact-match rate of each probe is tabulated by collection year,
to see whether the old probe loses matches as the binding site drifts while the new probe stays
stable. Subtype A records serve as the subtype-specificity control.

Output: h3_probe_audit_pilot_2026-10-28.tsv (probe x year x match rate)
"""
from collections import defaultdict, Counter
from pathlib import Path

import pandas as pd

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"

PROBES = {  # (name, 0-based start on the subtype B reference, forward sequence)
    "Old_N_5331890": ("AAGATGCAAATCATAAATTCACAGGA", 163),
    "New_N_5331890": ("GCATCATTGCTGTCATTAACTCAATT", 924),
}

# Alignment
aln, name = {}, None
for line in (WORK / "aln_N_full.fasta").read_text().split("\n"):
    if line.startswith(">"):
        name = line[1:].split()[0]
        aln[name] = []
    elif line.strip() and name is not None:
        aln[name].append(line.strip())
aln = {k: "".join(v).upper() for k, v in aln.items()}
ref_al = aln["REF_N"]
ref_cols = [i for i, c in enumerate(ref_al) if c != "-"]

# Year and subtype
mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["accession_base"] = mx["accession"].astype(str).str.split(".").str[0]
mx["year"] = mx["collection_date"].fillna("").str[-4:]
ok = mx["year"].str.match(r"^20\d\d$", na=False)
mx.loc[~ok, "year"] = mx["first_public"].fillna("").str[:4]
meta = mx.set_index("accession_base")["year"].to_dict()

st_of = {}
for line in (WORK / f"n_windows_full_{DAY}.fasta").read_text().split("\n"):
    if line.startswith(">"):
        p = line[1:].split()
        st_of[p[0]] = p[1].split("=")[1]

rows = []
for pname, (pseq, pos0) in PROBES.items():
    cols = ref_cols[pos0:pos0 + len(pseq)]
    stat = defaultdict(lambda: [0, 0])  # (st, year) -> [exact, n]
    statA = [0, 0]
    for accname, s in aln.items():
        if accname.startswith("REF_"):
            continue
        st = st_of.get(accname, "?")
        binding = "".join(s[i] for i in cols)
        if "N" in binding or "-" in binding:
            continue
        mm = sum(1 for a, b in zip(binding, pseq) if a != b)
        if st == "A":
            statA[0] += int(mm == 0)
            statA[1] += 1
            continue
        yr = meta.get(accname.split(".")[0], "")
        if not (yr and yr >= "2015" and yr <= "2026"):
            continue
        stat[(st, yr)][0] += int(mm == 0)
        stat[(st, yr)][1] += 1
    rows.append({"probe": pname, "subtype": "A", "year": "all_2015_2026",
                 "exact": statA[0], "n": statA[1],
                 "freq": round(statA[0] / max(1, statA[1]), 4)})
    for (st, yr), (e, n) in sorted(stat.items()):
        if n >= 30:
            rows.append({"probe": pname, "subtype": st, "year": yr,
                         "exact": e, "n": n, "freq": round(e / n, 4)})

df = pd.DataFrame(rows)
df.to_csv(WORK / f"h3_probe_audit_pilot_{DAY}.tsv", sep="\t", index=False)
print(df.to_string(index=False), flush=True)
