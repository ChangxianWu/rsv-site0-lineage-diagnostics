"""Position-resolved mismatch profile of the 17 audited arms (in silico stand-in for wet work).

For each arm, mismatches in the 2026 windows are binned by their position relative to the 3' end
of the oligonucleotide: -1 (terminal base), -2, -3, -4, -5, internal, and the 5'-terminal five
nucleotides.

The output reports frequencies only. It does not infer amplification impact, which requires the
wet-verification protocol and is not claimed here. Proportions are marginal, so one genome can
contribute to more than one position and the columns do not sum to the proportion with any
mismatch.

Output: h3_panel_position_profile_2026-10-28.tsv
"""
import time
from collections import Counter, defaultdict
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


def profile(label, arms, wins, subtype):
    rows = []
    for arm, seq in arms.items():
        n = 0
        pos = defaultdict(int)   # -1..-5, internal, 5p
        mm_total = 0
        for acc, w in wins.items():
            if year_of.get(acc.split(".")[0], "") != "2026":
                continue
            n += 1
            best = None
            for orient in ("fwd", "rev"):
                p_ = seq if orient == "fwd" else rc(seq)
                m = regex.search(f"(?:{p_}){{e<=2}}", w, flags=regex.BESTMATCH)
                if m is None:
                    continue
                seg = w[m.start():m.end()]
                mm = [i for i, (a, b) in enumerate(zip(seg, p_)) if a != b]
                if best is None or len(mm) < len(best[1]):
                    best = (orient, mm, len(p_))
                if not mm:
                    break
            if best is None:
                continue
            orient, mm, L = best
            mm_total += len(mm)
            for i in mm:
                d = (L - 1 - i) if orient == "fwd" else (L - 1 - (L - 1 - i))
                if orient == "rev":
                    d = i  # on a reverse-complement match the arm's 3' end is the segment start
                if 1 <= d <= 5:
                    pos[f"-{d}"] += 1
                elif i < 5 and orient == "fwd":
                    pos["5p"] += 1
                else:
                    pos["internal"] += 1
        if n == 0:
            continue
        row = {"gene": label, "oligo": arm, "subtype": subtype, "n_2026": n,
               "mismatch_events": mm_total}
        for k in ["-1", "-2", "-3", "-4", "-5", "internal", "5p"]:
            row[f"frac_pos_{k}"] = round(pos[k] / n, 4)
        rows.append(row)
    return rows


t0 = time.time()
rows = []
rows += profile("M", M_ARMS, read_windows(WORK / "win_M_full.fasta", "A"), "A")
rows += profile("L", L_ARMS, read_windows(WORK / "win_L_B.fasta"), "B")
rows += profile("N", N_PROBES, read_windows(WORK / "win_n_full.fasta", "B"), "B")
df = pd.DataFrame(rows)
df.to_csv(WORK / f"h3_panel_position_profile_{DAY}.tsv", sep="\t", index=False)
show = df.sort_values("frac_pos_-1", ascending=False)[
    ["gene", "oligo", "subtype", "n_2026", "frac_pos_-1", "frac_pos_-2", "frac_pos_-3",
     "frac_pos_internal", "frac_pos_5p"]]
print(show.to_string(index=False))
print("elapsed %.0fs" % (time.time() - t0))
