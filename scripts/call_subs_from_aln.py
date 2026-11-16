"""Call amino-acid substitutions in the reference coordinate frame.

Usage: python call_subs_from_aln.py <aln_fasta> <A|B> <out_tsv> <n_input_records>

The alignment produced by MAFFT --add --keeplength keeps the reference column set, so codons
are read only from reference non-gap columns and the frame cannot drift. Records with more than
30 gapped codons are counted and excluded from the denominator rather than silently scored.
"""
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

WORK = Path(__file__).resolve().parent
aln_file, st, out_tsv, n_in = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

F_COORDS = {"A": (5661, 7385), "B": (5666, 7390)}
ACC2ST = {"NC_038235.1": "A", "NC_001781.1": "B"}
COMP = str.maketrans("ACGTUN", "TGCAAN")

ref_f = {}
for block in (WORK / "refseq_rsv_fasta_2026-10-28.txt").read_text().split(">"):
    if not block.strip():
        continue
    acc = block.split("\n")[0].split()[0]
    if acc in ACC2ST:
        g = "".join(l.strip() for l in block.split("\n")[1:])
        a, b = F_COORDS[ACC2ST[acc]]
        ref_f[ACC2ST[acc]] = g[a - 1:b]

bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
i = 0
for x in bases:
    for y in bases:
        for z in bases:
            CODON[x + y + z] = aas[i]
            i += 1
REF_PROT = {st: "".join(CODON.get(ref_f[st][k:k + 3], "X")
                        for k in range(0, len(ref_f[st]) - 2, 3))[:-1] for st in ("A", "B")}

aln, name = {}, None
for line in (WORK / aln_file).read_text().split("\n"):
    if line.startswith(">"):
        name = line[1:].split()[0]
        aln[name] = []
    elif line.strip() and name is not None:
        aln[name].append(line.strip())
aln = {k: "".join(v).upper() for k, v in aln.items()}

ref_al = aln["REF_" + st]
ref_cols = [i for i, c in enumerate(ref_al) if c != "-"]
ref_nt = "".join(ref_al[i] for i in ref_cols)
assert ref_nt == ref_f[st], "ref frame mismatch"
codon_cols = [ref_cols[j:j + 3] for j in range(0, len(ref_cols) - 2, 3)]

subs = Counter()
n_rec = n_gap = 0
for accname, s in aln.items():
    if accname.startswith("REF_"):
        continue
    n_rec += 1
    gaps = 0
    aalist = []
    for cols in codon_cols:
        c = "".join(s[i] for i in cols)
        if "-" in c or "N" in c:
            gaps += 1
            aalist.append("-")
        else:
            aalist.append(CODON.get(c, "X"))
    if gaps > 30:
        n_gap += 1
        continue
    for j, aa in enumerate(aalist[:-1]):  # last group is the reference stop codon
        if aa != "-" and aa != REF_PROT[st][j]:
            subs[f"{REF_PROT[st][j]}{j + 1}{aa}"] += 1

den = n_rec - n_gap
rows = [(m, n, round(n / max(1, den), 5)) for m, n in subs.most_common()]
pd.DataFrame(rows, columns=["mutation", "n", "freq"]).to_csv(
    WORK / out_tsv, sep="\t", index=False, encoding="utf-8-sig")
print(f"aln={aln_file} subtype={st} input={n_in} parsed={n_rec} gap_heavy={n_gap} den={den}")
print("top15:")
for m, n, f in rows[:15]:
    print(" ", m, n, f)
