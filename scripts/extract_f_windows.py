"""Stage A: extraction of F-region windows by seed density.

For every complete genome of at least 13 kb, the reference F gene is indexed as 12-mers
(step 2). The (subtype, orientation) pair with the most hits is retained, the hit span is
taken as min/max of those hits, and the window is cut from 60 nt upstream to 72 nt
downstream of it. Windows are written with a subtype/orientation/hit-count header for
MAFFT --add --keeplength alignment.

A 2,000-nt sliding cluster selects the densest run of hits, so that scattered spurious
12-mer matches cannot inflate the span; windows whose span falls outside 1,500-1,950 nt
are rejected rather than aligned.

Outputs: f_windows_full_2026-10-28.fasta, f_windows_stats_2026-10-28.json
"""
import json
from collections import Counter
from pathlib import Path

WORK = Path(__file__).resolve().parent
RAW = WORK.parent / "raw"
FASTA = RAW / "h_rsv_nuccore_txid11250_2026-10-28.fasta"
DAY = "2026-10-28"
REFS = {"A": "NC_038235.1", "B": "NC_001781.1"}
F_COORDS = {"A": (5661, 7385), "B": (5666, 7390)}

COMP = str.maketrans("ACGTUN", "TGCAAN")
rc = lambda s: s.translate(COMP)[::-1]

# Reference F 12-mer index, both subtypes merged: key=(subtype, mer)
seed_map = {}
ref_f = {}
for st, acc in REFS.items():
    txt = (WORK / "refseq_rsv_fasta_2026-10-28.txt").read_text()
    for block in txt.split(">"):
        if block.startswith(acc):
            s = "".join(l.strip() for l in block.split("\n")[1:])
            a, b = F_COORDS[st]
            ref_f[st] = s[a - 1:b]
    for k in range(0, len(ref_f[st]) - 11, 2):
        seed_map.setdefault(ref_f[st][k:k + 12], []).append(st)
print("reference F lengths:", {k: len(v) for k, v in ref_f.items()},
      "seeds:", len(seed_map), flush=True)

records = []
acc, seq = None, []
with FASTA.open("r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith(">"):
            if acc is not None:
                records.append((acc, "".join(seq)))
            acc, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
records.append((acc, "".join(seq)))
print("records:", len(records), flush=True)

stats = Counter()
out = []
for acc, s in records:
    stats["total"] += 1
    if len(s) < 13000:
        continue
    t0 = s.upper().replace("U", "T")
    t0 = "".join(c if c in "ACGTN" else "N" for c in t0)
    t1 = rc(t0)
    # Hit positions per (subtype, orientation)
    hits = {"A": {0: [], 1: []}, "B": {0: [], 1: []}}
    for orient, t in ((0, t0), (1, t1)):
        for k in range(0, len(t) - 11):
            stl = seed_map.get(t[k:k + 12])
            if stl is not None:
                for st in stl:
                    hits[st][orient].append(k)
    best_st, best_orient, best_hits = None, None, []
    for st in ("A", "B"):
        for orient in (0, 1):
            if len(hits[st][orient]) > len(best_hits):
                best_st, best_orient, best_hits = st, orient, hits[st][orient]
    if len(best_hits) < 10:
        stats["seed_fail"] += 1
        continue
    stats[f"subtype_{best_st}"] += 1
    stats[f"orient_{'stored' if best_orient == 0 else 'rc'}"] += 1
    # Sliding-cluster density: densest run of hits within a 2,000-nt window
    hs = sorted(best_hits)
    best_i = 0
    best_cluster = (0, 0, 0)  # (count, lo, hi)
    j = 0
    for i in range(len(hs)):
        while hs[i] - hs[j] > 2000:
            j += 1
        if i - j + 1 > best_cluster[0]:
            best_cluster = (i - j + 1, hs[j], hs[i])
    cnt, lo, hi = best_cluster
    if cnt < 10 or not (1500 <= hi - lo <= 1950):
        stats["span_out"] += 1
        continue
    win_lo = max(0, lo - 60)
    win_hi = min(len(t0), hi + 72)
    window = (t0 if best_orient == 0 else t1)[win_lo:win_hi]
    span = hi - lo
    stats["extracted"] += 1
    out.append((acc, best_st, "stored" if best_orient == 0 else "rc",
                cnt, span, window))

with (WORK / f"f_windows_full_{DAY}.fasta").open("w") as fh:
    for acc, st, orient, nh, span, window in out:
        fh.write(f">{acc} subtype={st} orient={orient} seedhits={nh} span={span}\n")
        for k in range(0, len(window), 60):
            fh.write(window[k:k + 60] + "\n")

(WORK / f"f_windows_stats_{DAY}.json").write_text(
    json.dumps({k: int(v) for k, v in stats.items()}, indent=1))
print(json.dumps({k: int(v) for k, v in stats.items()}, indent=1), flush=True)
print("saved f_windows_full", flush=True)
