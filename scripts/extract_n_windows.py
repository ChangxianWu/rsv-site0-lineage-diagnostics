"""Stage A, N-gene version: extraction of N-region windows for the diagnostic pilot.

Reference N coding sequence: NC_038235 1140..2315 and NC_001781 1140..2315 (identical
coordinates, 1,176 nt each). The method matches the F-gene extractor: 12-mer seed density,
a 2,000-nt sliding cluster (1,400 nt here, because the N region is shorter), and a span
tolerance of 950-1,500 nt.

Outputs: n_windows_full_2026-10-28.fasta, n_windows_stats_2026-10-28.json
"""
import json
from collections import Counter
from pathlib import Path

WORK = Path(__file__).resolve().parent
RAW = WORK.parent / "raw"
FASTA = RAW / "h_rsv_nuccore_txid11250_2026-10-28.fasta"
DAY = "2026-10-28"
REFS = {"A": "NC_038235.1", "B": "NC_001781.1"}
N_COORDS = {"A": (1140, 2315), "B": (1140, 2315)}

COMP = str.maketrans("ACGTUN", "TGCAAN")
rc = lambda s: s.translate(COMP)[::-1]

seed_map = {}
ref_n = {}
for st, acc in REFS.items():
    for block in (WORK / "refseq_rsv_fasta_2026-10-28.txt").read_text().split(">"):
        if block.startswith(acc):
            s = "".join(l.strip() for l in block.split("\n")[1:])
            a, b = N_COORDS[st]
            ref_n[st] = s[a - 1:b]
    for k in range(len(ref_n[st]) - 11):
        seed_map.setdefault(ref_n[st][k:k + 12], []).append(st)
print("ref N lengths:", {k: len(v) for k, v in ref_n.items()},
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

stats = Counter()
out = []
for acc, s in records:
    stats["total"] += 1
    if len(s) < 13000:
        continue
    t0 = s.upper().replace("U", "T")
    t0 = "".join(c if c in "ACGTN" else "N" for c in t0)
    t1 = rc(t0)
    hits = {"A": {0: [], 1: []}, "B": {0: [], 1: []}}
    for orient, t in ((0, t0), (1, t1)):
        for k in range(len(t) - 11):
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
    hs = sorted(best_hits)
    best_i = j = 0
    best_cluster = (0, 0, 0)
    for i in range(len(hs)):
        while hs[i] - hs[j] > 1400:
            j += 1
        if i - j + 1 > best_cluster[0]:
            best_cluster = (i - j + 1, hs[j], hs[i])
    cnt, lo, hi = best_cluster
    if cnt < 10 or not (950 <= hi - lo <= 1500):
        stats["span_out"] += 1
        continue
    window = (t0 if best_orient == 0 else t1)[max(0, lo - 60):min(len(t0), hi + 72)]
    stats["extracted"] += 1
    out.append((acc, best_st, "stored" if best_orient == 0 else "rc", cnt, hi - lo, window))

with (WORK / f"n_windows_full_{DAY}.fasta").open("w") as fh:
    for acc, st, orient, nh, span, window in out:
        fh.write(f">{acc} subtype={st} orient={orient} seedhits={nh} span={span}\n")
        for k in range(0, len(window), 60):
            fh.write(window[k:k + 60] + "\n")
(WORK / f"n_windows_stats_{DAY}.json").write_text(
    json.dumps({k: int(v) for k, v in stats.items()}, indent=1))
print(json.dumps({k: int(v) for k, v in stats.items()}, indent=1), flush=True)
