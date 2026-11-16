"""Reference-free audit of published diagnostic oligonucleotides, and Figure 4.

Method: each arm is searched as an exact string in the orientation-normalised gene windows
(a forward or reverse-complement hit counts), tabulated by collection year. Searching the window
itself rather than a reference-aligned frame removes the reference-era bias that made
contemporary arms look absent from a 1956-frame alignment. The N-gene probes reuse the
frame-based result, which matches the subtype B reference exactly and is therefore equivalent.

Classification uses the absolute 2026 match rate (< 80% deep mismatch, 80-90% declining,
>= 90% stable) so that the category cannot be steered by the choice of baseline years.

Outputs: h3_panel_audit_unified_2026-10-28.tsv, h3_panel_audit_summary_unified_2026-10-28.tsv
and Figure4.(pdf|tif) for the 17 audited arms.
"""
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

WORK = Path(__file__).resolve().parent
OUT = WORK / "figures"
OUT.mkdir(exist_ok=True)
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

mx = pd.read_csv(WORK / "ena_rsv_meta_matrix_sp11250_2026-10-27.tsv", sep="\t", dtype=str)
mx["acc"] = mx.accession.astype(str).str.split(".").str[0]
yr = mx.collection_date.fillna("").str[-4:]
ok = yr.str.match(r"^20\d\d$", na=False)
mx["year"] = yr.where(ok, mx.first_public.fillna("").str[:4])
year_of = mx.set_index("acc")["year"].to_dict()


def read_windows(path, only_subtype=None):
    wins, name, buf = {}, None, []
    keep = True
    for line in Path(path).read_text().split("\n"):
        if line.startswith(">"):
            if name and keep:
                wins[name] = "".join(buf)
            p = line[1:].split()
            name, buf = p[0], []
            keep = True
            if only_subtype:
                if "=" in line:
                    tags = dict(kv.split("=") for kv in p[1:] if "=" in kv)
                    keep = tags.get("subtype") == only_subtype
                else:  # header like ">ACC A"
                    keep = len(p) > 1 and p[1] == only_subtype
        elif line.strip():
            buf.append(line.strip())
    if name and keep:
        wins[name] = "".join(buf)
    return wins


def audit(label, arms, wins):
    rows = []
    for arm, seq in arms.items():
        variants = [seq, rc(seq)]
        stat = defaultdict(lambda: [0, 0])
        for acc, w in wins.items():
            y = year_of.get(acc.split(".")[0], "")
            if not (y and "2015" <= y <= "2026"):
                continue
            hit = any(v in w for v in variants)
            stat[y][0] += int(hit)
            stat[y][1] += 1
        for y, (h, n) in sorted(stat.items()):
            if n >= 30:
                rows.append({"gene": label, "arm": arm, "year": y,
                             "hits": h, "n": n, "rate": h / n})
    return rows


# M windows are orientation-normalised; subtype A records only (the M assay targets subtype A)
wins_m = read_windows(WORK / "win_M_full.fasta", only_subtype="A")
rows = audit("M", M_ARMS, wins_m)
wins_l = read_windows(WORK / "win_L_B.fasta")
rows += audit("L", L_ARMS, wins_l)

df = pd.DataFrame(rows)
df.to_csv(WORK / f"h3_panel_audit_unified_{DAY}.tsv", sep="\t", index=False)

# N-gene probes: reuse the frame-based result
np_ = pd.read_csv(WORK / f"h3_probe_audit_pilot_{DAY}.tsv", sep="\t")
np_ = np_[(np_.subtype == "B") & np_.year.astype(str).str.match(r"^\d{4}$")]
np_["year"] = np_.year.astype(int).astype(str)
np_["gene"] = "N"; np_["arm"] = np_.probe
np_ = np_.rename(columns={"exact": "hits"})[["gene", "arm", "year", "hits", "n", "freq"]]
np_ = np_.rename(columns={"freq": "rate"})
np_["arm"] = np_["arm"].str.replace("byyear_B", "", regex=False).str.rstrip("_")

all_df = pd.concat([df[["gene", "arm", "year", "hits", "n", "rate"]],
                    np_[["gene", "arm", "year", "hits", "n", "rate"]]], ignore_index=True)

summary = []
for (g, a), d in all_df.groupby(["gene", "arm"]):
    y26 = d[d.year == "2026"]
    base = d[d.year.isin(["2015", "2016", "2017", "2018"])]
    if y26.n.sum() < 30:
        continue
    r26 = y26.hits.sum() / y26.n.sum()
    rb = base.hits.sum() / base.n.sum() if base.n.sum() >= 30 else np.nan
    delta = r26 - rb
    shape = ("deep mismatch" if r26 < 0.80 else
             "declining" if r26 < 0.90 else "stable")
    summary.append({"gene": g, "oligo": a, "rate_2026": round(r26, 3),
                    "baseline_2015_18": (round(rb, 3) if not np.isnan(rb) else None),
                    "delta": (round(delta, 3) if not np.isnan(delta) else None),
                    "signal_shape": shape,
                    "subtype": "B" if g in ("N", "L") else "A"})
S = pd.DataFrame(summary).sort_values("rate_2026").reset_index(drop=True)
S.to_csv(WORK / f"h3_panel_audit_summary_unified_{DAY}.tsv", sep="\t", index=False)
print(S.to_string(index=False))

# Figure 4
marker = {"stable": "o", "declining": "s", "deep mismatch": "^"}
fig, ax = plt.subplots(figsize=(7.0, 5.4))
y = np.arange(len(S))
for i, r in S.iterrows():
    ax.plot(r.rate_2026 * 100, i, marker[r.signal_shape], color="black", ms=5.5,
            mfc="white" if r.signal_shape == "stable" else "black")
    ax.text(r.rate_2026 * 100 + 1.3, i, f"{r.rate_2026*100:.1f}%", va="center", fontsize=6,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
ax.axvline(95, color="black", ls=":", lw=0.8)
ax.set_yticks(y, [f"{r.oligo} ({r.gene}, {r.subtype})" for _, r in S.iterrows()], fontsize=6.4)
ax.set_xlabel("Exact-match rate to 2026 circulating genomes (%)")
ax.set_xlim(0, 112)
ax.set_title("Published diagnostic oligonucleotides versus 2026 RSV genomes", fontsize=9)
from matplotlib.lines import Line2D
leg = [Line2D([0], [0], marker=marker[k], color="black", ls="", ms=5,
              mfc="white" if k == "stable" else "black", label=k)
       for k in ["deep mismatch", "declining", "stable"]]
ax.legend(handles=leg, fontsize=6.3, loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig(OUT / "Figure4.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "Figure4.tif", dpi=300, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print("Figure4 regenerated (17 arms)")
print("range 2026: %.1f%% - %.1f%%" % (S.rate_2026.min() * 100, S.rate_2026.max() * 100))
