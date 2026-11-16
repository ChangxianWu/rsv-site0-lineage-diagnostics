"""H2 trajectory layer: key-substitution frequencies by year and country group.

Inputs: aln_A_full.fasta / aln_B_full.fasta (reference coordinate frame from MAFFT
--keeplength) and the ENA metadata matrix (country, collection_date, first_public).

Method: substitutions are called per record within the site-O window (positions 1-260 for both
subtypes), joined to metadata by accession base, dated by collection year with first_public year
as fallback, and aggregated. Country groups compare the ten anchor exposure countries
(ESP, FRA, USA, LUX, AUT, IRL, DEU, NLD, BEL, PRT) against other classifiable countries.

Outputs: traj_key_sites_by_year_2026-10-28.tsv (mutation x year x group x n x den x freq) and
traj_per_record_site0_2026-10-28.tsv.gz (sparse per-record substitutions, reused by the lineage
and diagnostic layers).
"""
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
F_COORDS = {"A": (5661, 7385), "B": (5666, 7390)}
ACC2ST = {"NC_038235.1": "A", "NC_001781.1": "B"}
ANCHORS = {"ESP", "FRA", "USA", "LUX", "AUT", "IRL", "DEU", "NLD", "BEL", "PRT"}
KEYMUTS = {"A": {"K65R", "E66K", "N88S", "A74T", "K68E", "K68Q", "K68N"},
           "B": {"K68Q", "K68N", "I206M", "Q209R", "S211N", "S190N", "K191R"}}
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

bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
i = 0
for x in bases:
    for y in bases:
        for z in bases:
            CODON[x + y + z] = aas[i]
            i += 1

ref_f = {}
for block in (WORK / "refseq_rsv_fasta_2026-10-28.txt").read_text().split(">"):
    if not block.strip():
        continue
    acc = block.split("\n")[0].split()[0]
    if acc in ACC2ST:
        g = "".join(l.strip() for l in block.split("\n")[1:])
        a, b = F_COORDS[ACC2ST[acc]]
        ref_f[ACC2ST[acc]] = g[a - 1:b]
REF_PROT = {st: "".join(CODON.get(ref_f[st][k:k + 3], "X")
                        for k in range(0, len(ref_f[st]) - 2, 3))[:-1] for st in ("A", "B")}

# Metadata
mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["accession_base"] = mx["accession"].astype(str).str.split(".").str[0]
mx["iso3"] = mx["country"].fillna("").str.split(":").str[0].str.strip().map(CN).fillna("")
mx["year"] = mx["collection_date"].fillna("").str[-4:]
mx.loc[~mx["year"].str.match(r"^20\d\d$", na=False), "year"] = \
    mx["first_public"].fillna("").str[:4]
meta = mx.set_index("accession_base")[["iso3", "year"]].to_dict("index")

# Sparse per-record substitutions (positions 1-260) plus key-mutation aggregation
rows_out = []
agg = defaultdict(Counter)      # (st, mut, year, group) -> n
den = defaultdict(Counter)      # (st, year, group) -> den
gz = gzip.open(WORK / f"traj_per_record_site0_{DAY}.tsv.gz", "wt", encoding="utf-8")
gz.write("accession\tsubtype\tyear\tiso3\tgroup\tmutations\n")

for st in ("A", "B"):
    aln, name = {}, None
    for line in (WORK / f"aln_{st}_full.fasta").read_text().split("\n"):
        if line.startswith(">"):
            name = line[1:].split()[0]
            aln[name] = []
        elif line.strip() and name is not None:
            aln[name].append(line.strip())
    aln = {k: "".join(v).upper() for k, v in aln.items()}
    ref_al = aln["REF_" + st]
    ref_cols = [i for i, c in enumerate(ref_al) if c != "-"]
    codon_cols = [ref_cols[j:j + 3] for j in range(0, len(ref_cols) - 2, 3)][:len(REF_PROT[st])]
    for accname, s in aln.items():
        if accname.startswith("REF_"):
            continue
        aalist = []
        for cols in codon_cols:
            c = "".join(s[i] for i in cols)
            aalist.append("-" if ("-" in c or "N" in c) else CODON.get(c, "X"))
        if aalist.count("-") > 30:
            continue
        m = meta.get(accname.split(".")[0], {})
        year = m.get("year", "")
        iso3 = m.get("iso3", "")
        if not year or year < "2015" or year > "2026":
            continue
        group = "anchor" if iso3 in ANCHORS else ("other" if iso3 else "na")
        if group == "na":
            continue
        den[(st, year, group)]["den"] += 1
        muts = []
        for j, aa in enumerate(aalist):
            if aa != "-" and aa != REF_PROT[st][j]:
                mut = f"{REF_PROT[st][j]}{j + 1}{aa}"
                muts.append(mut)
                if j < 260:
                    agg[(st, mut, year, group)]["n"] += 1
        gz.write(f"{accname}\t{st}\t{year}\t{iso3}\t{group}\t{';'.join(muts)}\n")
gz.close()

rows = []
for (st, mut, year, group), c in agg.items():
    d = den[(st, year, group)]["den"]
    if mut in KEYMUTS[st] and d >= 30:
        rows.append({"subtype": st, "mutation": mut, "year": year, "group": group,
                     "n": c["n"], "den": d, "freq": round(c["n"] / d, 4)})
df = pd.DataFrame(rows).sort_values(["subtype", "mutation", "group", "year"])
df.to_csv(WORK / f"traj_key_sites_by_year_{DAY}.tsv", sep="\t", index=False, encoding="utf-8-sig")
print(df.to_string(index=False), flush=True)
print("saved traj_key_sites_by_year + per_record gz", flush=True)
