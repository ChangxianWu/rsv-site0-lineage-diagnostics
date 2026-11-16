"""Mismatch-distance distribution and 3'-proximal mismatch rates per oligonucleotide arm.

Exact-match rate is not detection performance. For every arm the window is searched with a fuzzy
matcher allowing up to two errors in either orientation, and four quantities are tabulated:
p_exact, p_le1, p_le2, p_ge3, plus the proportion carrying at least one mismatch within the
3'-terminal five nucleotides of the arm.

Both the 2026 window and the 2015-2018 baseline are computed, so drift can be read as a change in
distance distribution rather than as a single rate.

Output: h3_panel_mismatch_dist_2026-10-28.tsv
"""
import time
from pathlib import Path

import pandas as pd
import regex

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
COMP = str.maketrans("ACGTUN", "TGCAAN")
rc = lambda s: s.translate(COMP)[::-1]

M_ARMS = {"F3": "GCTGTTCAATACAATGTCCTAGA", "B3": "GGTAAATTTGCTGGGCATT",
          "FIP_F1c": "TCTGCTGGCATGGATGATTGG", "FIP_F2": "AGACGATGATCCTGCATCA",
          "BIP_B1c": "CTAGTGAAACAAATATCCACACCCA", "BIP_B2": "GCACTGCACTTCTTGAGTT",
          "LF": "ACATGGGCACCCATATTGTAAG"}
L_ARMS = {"F3": "TTTGGTGGTGGTGATCCT", "B3": "CAAACTCGGCATTTGGATT",
          "FIP_F1c": "AGCTCAACACAAACACTGAATGTAC", "FIP_F2": "TATATCGAAGCTTTTATAGGAGAAC",
          "BIP_B1c": "CTGGTCACGATCTACAAGATAAGC", "BIP_B2": "TGTGATGACACAAGTCAAGAA",
          "LF": "TATAGCTTCTGTAAGGAAGTCTGGA", "LB": "GGATCTTCCAGATGATAGACTGAAC"}
# N-gene probes (searching subtype B N windows; a probe is not a primer, so its 5'/3' ends are
# treated as the probe's own termini)
N_PROBES = {"Old_N_5331890": "AAGATGCAAATCATAAATTCACAGGA",
            "New_N_5331890": "GCATCATTGCTGTCATTAACTCAATT"}

mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["acc"] = mx.accession.astype(str).str.split(".").str[0]
yr = mx.collection_date.fillna("").str[-4:]
ok = yr.str.match(r"^20\d\d$", na=False)
mx["year"] = yr.where(ok, mx.first_public.fillna("").str[:4])
year_of = mx.set_index("acc")["year"].to_dict()


def read_windows(path, only_subtype=None):
    wins, name, buf, keep = {}, None, [], True
    for line in Path(path).read_text().split("\n"):
        if line.startswith(">"):
            if name and keep:
                wins[name] = "".join(buf)
            p = line[1:].split()
            name, buf, keep = p[0], [], True
            if only_subtype:
                if any("=" in x for x in p[1:]):
                    tags = dict(kv.split("=") for kv in p[1:] if "=" in kv)
                    keep = tags.get("subtype") == only_subtype
                else:
                    keep = len(p) > 1 and p[1] == only_subtype
        elif line.strip():
            buf.append(line.strip())
    if name and keep:
        wins[name] = "".join(buf)
    return wins


def best_match(w, pat):
    """Return (distance, mismatch positions counted from the 5' end of pat), both orientations."""
    best = (99, [])
    for orient in ("fwd", "rev"):
        p = pat if orient == "fwd" else rc(pat)
        m = regex.search(f"(?:{p}){{e<=2}}", w, flags=regex.BESTMATCH)
        if m is None:
            continue
        seg = w[m.start():m.end()]
        mm = [i for i, (a, b) in enumerate(zip(seg, p)) if a != b]
        if len(mm) < best[0]:
            # On a reverse-complement match the 3' end of pat is the 5' end of p
            best = (len(mm), mm if orient == "fwd" else [len(p) - 1 - i for i in mm])
        if best[0] == 0:
            break
    return best


def audit(label, arms, wins, subtype, years):
    rows = []
    for arm, seq in arms.items():
        accs = [(a, w) for a, w in wins.items()
                if year_of.get(a.split(".")[0], "") in years]
        d0 = d1 = d2 = d3 = prox = 0
        for a, w in accs:
            dist, mmpos = best_match(w, seq)
            d0 += dist == 0
            d1 += dist <= 1
            d2 += dist <= 2
            d3 += dist >= 3
            if dist <= 2 and any(p >= len(seq) - 5 for p in mmpos):
                prox += 1
        n = len(accs)
        if n == 0:
            continue
        rows.append({"gene": label, "oligo": arm, "subtype": subtype,
                     "window": "2026" if years == ("2026",) else "2015-18",
                     "n": n, "p_exact": round(d0 / n, 4), "p_le1": round(d1 / n, 4),
                     "p_le2": round(d2 / n, 4), "p_ge3": round(d3 / n, 4),
                     "p_3prox": round(prox / n, 4)})
    return rows


t0 = time.time()
rows = []
rows += audit("M", M_ARMS, read_windows(WORK / "win_M_full.fasta", only_subtype="A"), "A", ("2026",))
rows += audit("L", L_ARMS, read_windows(WORK / "win_L_B.fasta"), "B", ("2026",))
rows += audit("N", N_PROBES, read_windows(WORK / "win_n_full.fasta"), "B", ("2026",))
# Baseline (2015-18)
rows += audit("M", M_ARMS, read_windows(WORK / "win_M_full.fasta", only_subtype="A"), "A",
              ("2015", "2016", "2017", "2018"))
rows += audit("L", L_ARMS, read_windows(WORK / "win_L_B.fasta"), "B",
              ("2015", "2016", "2017", "2018"))
rows += audit("N", N_PROBES, read_windows(WORK / "win_n_full.fasta"), "B",
              ("2015", "2016", "2017", "2018"))

df = pd.DataFrame(rows)
df.to_csv(WORK / f"h3_panel_mismatch_dist_{DAY}.tsv", sep="\t", index=False)
p = df[df.window == "2026"].sort_values("p_exact")
print(p[["gene", "oligo", "subtype", "n", "p_exact", "p_le1", "p_le2", "p_ge3", "p_3prox"]].to_string(index=False))
print(f"elapsed {time.time()-t0:.0f}s")
