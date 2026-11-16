"""Phylogenetic check on monophyly.

Step "prep": stratified subsample of subtype B records (450 B.D.E.1-marker-consistent, 250
B.D.4.1.1-only, 300 other, selected deterministically by sha256 of the accession, plus the
subtype B reference) written as a FASTA subset together with a trait table for the six lineage
markers.

Step "fitch": parse the IQ-TREE tree and count independent origins of each marker under Fitch
parsimony (post-order set intersection, then a pre-order pass that charges one change whenever
the parent state is absent from a child's set). The root state is taken by majority; no
per-branch bootstrap threshold is applied.

Usage: python phylo_ed3.py prep | python phylo_ed3.py fitch
"""
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

WORK = Path(__file__).resolve().parent
DAY = "2026-10-28"

# Lineage marker sets
BD411 = {"K191R", "I206M", "Q209R"}
BDE1 = {"S190N", "S211N", "S389P"}


def prep():
    # Read the subtype B alignment (F region, 1,725 nt reference frame)
    seqs, name, buf = {}, None, []
    for line in (WORK / "aln_B_full.fasta").read_text().split("\n"):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()

    # Per-record substitutions (from the sparse per-record table)
    import pandas as pd
    per = pd.read_csv(WORK / f"traj_per_record_site0_{DAY}.tsv.gz", sep="\t", dtype=str)
    per = per[per.subtype == "B"]
    muts = dict(zip(per.accession, per.mutations.fillna("")))

    groups = {"bde1": [], "bd411_only": [], "other": []}
    for acc, s in seqs.items():
        if acc.startswith("REF_"):
            continue
        m = set(muts.get(acc, "").split(";")) - {""}
        if not m:
            continue
        if BD411 <= m and BDE1 <= m:
            groups["bde1"].append(acc)
        elif BD411 <= m:
            groups["bd411_only"].append(acc)
        else:
            groups["other"].append(acc)
    print({k: len(v) for k, v in groups.items()})

    n_each = {"bde1": 450, "bd411_only": 250, "other": 300}
    picked = []
    for g, n in n_each.items():
        pool = sorted(groups[g], key=lambda a: hashlib.sha256(a.encode()).hexdigest())
        picked += pool[:n]
    picked.append("REF_B")
    with (WORK / "tree_B_subset.fasta").open("w") as fh:
        for a in picked:
            s = seqs[a]
            fh.write(f">{a}\n")
            for k in range(0, len(s), 60):
                fh.write(s[k:k + 60] + "\n")

    # Trait table (the F region spans residues 1-574, so S389P is inside the window)
    trait = {}
    for a in picked:
        if a == "REF_B":
            trait[a] = {k: 0 for k in ["K191R", "I206M", "Q209R", "S190N", "S211N", "S389P"]}
            continue
        m = set(muts.get(a, "").split(";")) - {""}
        trait[a] = {k: int(k in m) for k in
                    ["K191R", "I206M", "Q209R", "S190N", "S211N", "S389P"]}
    (WORK / "tree_B_traits.json").write_text(json.dumps(trait, indent=0))
    print("subset n =", len(picked))


def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        children = []
        name = ""
        if s[pos] == "(":
            pos += 1
            while True:
                children.append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        m = re.match(r"[^,():;]+", s[pos:])
        if m:
            name = m.group(0)
            pos += len(name)
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[0-9.eE\-+]+", s[pos:])
            if m:
                pos += len(m.group(0))
        return {"name": name.strip(), "children": children}

    return node()


def fitch_origins(tree, states):
    """Return the number of changes for one binary trait (Fitch post-order + pre-order count)."""
    def post(n):
        if not n["children"]:
            st = states.get(n["name"])
            n["_set"] = {st} if st in (0, 1) else {0, 1}
            return n["_set"]
        sets = [post(c) for c in n["children"]]
        inter = set.intersection(*sets)
        n["_set"] = inter if inter else set.union(*sets)
        return n["_set"]

    def count(n, parent_state=None):
        s = n["_set"]
        if parent_state is None:
            state = next(iter(s)) if len(s) == 1 else 0
        else:
            state = parent_state if parent_state in s else next(iter(s))
        total = 0
        if n["children"]:
            for c in n["children"]:
                cs = c["_set"]
                if state in cs:
                    total += count(c, state)
                else:
                    total += 1 + count(c, next(iter(cs)))
        n["_state"] = state
        return total

    post(tree)
    return count(tree)


def fitch():
    trait = json.loads((WORK / "tree_B_traits.json").read_text())
    tree_line = ""
    for line in (WORK / "tree_B.treefile").read_text().split("\n"):
        if line.strip().startswith("("):
            tree_line = line.strip()
            break
    if not tree_line:
        print("no tree")
        return
    rows = []
    for marker in ["K191R", "I206M", "Q209R", "S190N", "S211N", "S389P"]:
        states = {a: v[marker] for a, v in trait.items()}
        changes = fitch_origins(parse_newick(tree_line), states)
        n1 = sum(states.values())
        rows.append({"trait": marker, "n_carrying": n1, "independent_origins_parsimony": changes})
    # Combination traits: carriage of all three markers of each set
    states = {a: int(all(v[k] == 1 for k in BDE1)) for a, v in trait.items()}
    changes = fitch_origins(parse_newick(tree_line), states)
    rows.append({"trait": "BDE1_3marker_combo", "n_carrying": sum(states.values()),
                 "independent_origins_parsimony": changes})
    states = {a: int(all(v[k] == 1 for k in BD411)) for a, v in trait.items()}
    changes = fitch_origins(parse_newick(tree_line), states)
    rows.append({"trait": "BD411_3marker_combo", "n_carrying": sum(states.values()),
                 "independent_origins_parsimony": changes})
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(WORK / f"phylo_marker_origins_{DAY}.tsv", sep="\t", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    {"prep": prep, "fitch": fitch}[sys.argv[1]]()
