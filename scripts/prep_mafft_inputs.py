"""Prepare MAFFT inputs (pure file operations, no subprocess).

Outputs: ref_A.fasta / ref_B.fasta (reference F nucleotides, the coordinate anchor),
win_A_pilot.fasta (first 1,000 subtype A windows), win_A_full.fasta and win_B_full.fasta.
Windows flagged orient=rc are reverse-complemented to the stored reference orientation.
"""
from pathlib import Path

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"
F_COORDS = {"A": (5661, 7385), "B": (5666, 7390)}  # NC_038235 / NC_001781
ACC2ST = {"NC_038235.1": "A", "NC_001781.1": "B"}
COMP = str.maketrans("ACGTUN", "TGCAAN")
rc = lambda s: s.translate(COMP)[::-1]

ref_f = {}
for block in (WORK / "refseq_rsv_fasta_2026-10-28.txt").read_text().split(">"):
    if not block.strip():
        continue
    acc = block.split("\n")[0].split()[0]
    if acc in ACC2ST:
        st = ACC2ST[acc]
        g = "".join(l.strip() for l in block.split("\n")[1:])
        a, b = F_COORDS[st]
        ref_f[st] = g[a - 1:b]
        (WORK / f"ref_{st}.fasta").write_text(f">REF_{st}\n{ref_f[st]}\n")

wins = {"A": [], "B": []}
cur = None
buf = []
for line in (WORK / f"f_windows_full_{DAY}.fasta").read_text().split("\n"):
    if line.startswith(">"):
        if cur is not None:
            wins[cur[1]].append((cur[0], cur[2], "".join(buf)))
        p = dict(kv.split("=") for kv in line[1:].split()[1:])
        cur = (line[1:].split()[0], p["subtype"], p["orient"])
        buf = []
    elif line.strip() and cur is not None:
        buf.append(line.strip())
if cur is not None:
    wins[cur[1]].append((cur[0], cur[2], "".join(buf)))

def dump(name, subset):
    with (WORK / name).open("w") as fh:
        for acc, orient, s in subset:
            seq = s if orient == "stored" else rc(s)
            fh.write(f">{acc}\n")
            for k in range(0, len(seq), 60):
                fh.write(seq[k:k + 60] + "\n")

dump("win_A_pilot.fasta", wins["A"][:1000])
dump("win_A_full.fasta", wins["A"])
dump("win_B_full.fasta", wins["B"])
print("A total:", len(wins["A"]), "| B total:", len(wins["B"]),
      "| pilot:", min(1000, len(wins["A"])), flush=True)
print("ref F lengths:", {k: len(v) for k, v in ref_f.items()}, flush=True)
