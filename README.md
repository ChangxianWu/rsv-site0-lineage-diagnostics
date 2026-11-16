# Lineage replacement and diagnostic oligonucleotide divergence in RSV

Analysis code for the manuscript:

> **Lineage Replacement Underlies Site Ø Substitution Trends in Respiratory Syncytial Virus and
> Diagnostic Oligonucleotide Divergence, 23 Countries, 2015–2026**

Submitted to the *Journal of Medical Virology*. Sole author and corresponding author:
Changxian Wu, College of Agriculture and Animal Husbandry Science and Technology, Ningxia
Vocational Technical College for Minorities, Wuzhong 751100, Ningxia, China
(wchangx173@163.com).

This repository contains the analysis code only: no sequence data, no derived tables, no
credentials. Every sequence analysed is public and is re-downloaded by the scripts below.

---

## What the code does

| Layer | Question | Script |
|---|---|---|
| Retrieval | all public human RSV nucleotide records | `download_ed3_nuccore_fasta.py` |
| Windows | F-region windows by 12-mer seed density | `extract_f_windows.py` |
| Windows | N-region windows (diagnostic probes) | `extract_n_windows.py` |
| Alignment input | reference frame + orientation normalisation | `prep_mafft_inputs.py` |
| Variant calling | codon-aware substitutions in the reference frame | `call_subs_from_aln.py` |
| Trajectories | substitution frequency by year and country group | `traj_site0_by_year.py` |
| Seasons | collection season and the analysis base | `seasonize_and_h1base.py`, `build_analysis_base_v1.py` |
| Exposure model | primary model, sensitivity battery, bootstrap, permutation | `h1_frozen_outcome_run.py`, `h1_sensitivity_battery.py` |
| Phylogeny | independent origins of lineage markers (Fitch parsimony) | `phylo_ed3.py` |
| Diagnostics | exact/≤1/≤2 mismatch, 3′-proximal position profiles, Figure 4 | `h3_panel_audit_unified.py`, `h3_panel_mismatch_dist.py`, `h3_panel_position_profile.py`, `h3_probe_audit_pilot.py` |

Two things this code deliberately does **not** do: it generates no neutralisation or
amplification phenotypes (all diagnostic statements in the paper are in silico sequence
compatibility), and it does not use GISAID data.

## Environment

```
python >= 3.10
pip install -r requirements.txt
mafft >= 7.5      # --add --keeplength
iqtree2 >= 2.2    # GTR+F+R2, --ufboot 1000 --bnni
```

`requirements.txt` covers requests, pandas, numpy, statsmodels, regex and matplotlib.
MAFFT and IQ-TREE are external binaries and are not pip-installable.

## Inputs

**NCBI Nucleotide (fetched by script).** `download_ed3_nuccore_fasta.py` queries
`txid11250[Organism]` through E-utilities, batches efetch at 200 UIDs, resumes between batches and
writes a manifest with the file sha256. The taxonomy identifier is resolved from the scientific
name on the day of retrieval: the same numeric identifier has been confused with rabies virus in
this project's history, so always confirm that the returned `scientific_name` values contain
"respiratory syncytial" before analysing anything.

**ENA metadata matrix (not redistributed, 9.7 MB).** One call to the ENA Portal API reproduces it:

```
https://www.ebi.ac.uk/ena/portal/api/search
  ?result=sequence
  &query=tax_tree(12814)
  &fields=accession,tax_id,scientific_name,country,collection_date,host,isolate,strain,
          keywords,mol_type,sample_accession,study_accession,first_public,base_count
  &format=tsv&limit=0
```

Save as `scripts/ena_rsv_meta_matrix_sp11250_2026-10-27.tsv` (tab-separated, header included).
The matrix carries 71,842 rows; country is recorded for 94.9% and collection date for 90.6%.

**Reference sequences.** `scripts/refseq_rsv_fasta_2026-10-28.txt` holds the two RefSeq genomes
that anchor the coordinate frame: NC_038235.1 (subtype A) and NC_001781.1 (subtype B).

## Running the pipeline

All scripts resolve their working directory from `Path(__file__)`, so run them from `scripts/`
and they read and write there; the downloaded FASTA goes to `../raw/`.

```bash
cd scripts

# 1. sequences (~74k records, resumable, NCBI rate limit respected)
python download_ed3_nuccore_fasta.py

# 2. windows
python extract_f_windows.py        # -> f_windows_full_2026-10-28.fasta
python extract_n_windows.py        # -> n_windows_full_2026-10-28.fasta

# 3. reference frame and orientation
python prep_mafft_inputs.py        # -> ref_A.fasta, ref_B.fasta, win_A_full.fasta, win_B_full.fasta

# 4. reference-frame alignment (one per subtype; N region likewise against ref_N)
mafft --add win_A_full.fasta --keeplength --anysymbol --thread -1 ref_A.fasta > aln_A_full.fasta
mafft --add win_B_full.fasta --keeplength --anysymbol --thread -1 ref_B.fasta > aln_B_full.fasta

# 5. substitutions in the reference frame
python call_subs_from_aln.py aln_A_full.fasta A subs_A_full_2026-10-28.tsv 13320
python call_subs_from_aln.py aln_B_full.fasta B subs_B_full_2026-10-28.tsv 11643

# 6. trajectories, seasons, analysis base
python traj_site0_by_year.py
python seasonize_and_h1base.py
python build_analysis_base_v1.py

# 7. exposure models
python h1_frozen_outcome_run.py
python h1_sensitivity_battery.py

# 8. phylogeny (F+G concatenation at 2,625 nt, then IQ-TREE with 1,000 UFBoot)
python phylo_ed3.py prep
iqtree2 -s tree_B_subset.fasta -m GTR+F+R2 --ufboot 1000 --bnni -T 4 --prefix tree_B
python phylo_ed3.py fitch

# 9. diagnostic oligonucleotide audit
python h3_panel_audit_unified.py
python h3_panel_mismatch_dist.py
python h3_panel_position_profile.py
python h3_probe_audit_pilot.py
```

### File-name map

The assay layers use their own working names, so a few intermediates are renamed between steps.
Nothing is recomputed; these are the same files under the names each script expects:

| Produced by | Renamed to | Used as |
|---|---|---|
| `prep_mafft_inputs.py` → `win_A_full.fasta` | `win_M_full.fasta` | subtype A multiplex arms (M) |
| `prep_mafft_inputs.py` → `win_B_full.fasta` | `win_L_B.fasta` | subtype B LAMP arms (L) |
| `extract_n_windows.py` → `n_windows_full_2026-10-28.fasta` | `win_n_full.fasta` | N-gene probes (N) |
| MAFFT output for the N region | `aln_N_full.fasta` (reference row named `REF_N`) | probe audit |

The alignment scripts expect the reference record to be named `REF_A`, `REF_B` or `REF_N`; rename
the row after MAFFT, which carries over the name from `ref_X.fasta`.

## Reproducibility notes

- **Determinism.** Subsampling is ordered by `sha256(accession)` rather than by an RNG; the
  bootstrap and permutation use a fixed seed (`20261028`). Re-running produces the same tables.
- **Denominators.** Marker frequencies use subtype B records with interpretable codons at every
  marker position and a collection year in 2015–2026 (n = 10,174). Escape-candidate frequencies
  use the aligned comparison set (n = 11,643). Subtype A comparisons use 13,320 records. The
  marker sets are nested within the 10,174 denominator: B.D.E.1 6,269 ⊆ B.D.4.1.1 8,665 ⊆ I206M
  carriers 9,127.
- **Gap handling.** Records with more than 30 gapped codons are excluded from the denominator
  rather than scored as wild type; the count is printed by `call_subs_from_aln.py`.
- **Reference-era bias.** Diagnostic arms are searched directly in the orientation-normalised
  window, not against a reference-aligned frame. Against the 1956 subtype B reference, three arms
  falsely appear to have no exact match; the window search restores 92–96% compatibility.
- **Aggregate tables.** The tables that back every reported value are the Supplementary Tables
  S1–S10 of the manuscript; this repository regenerates them from public data.

## Changes made for release

The scripts are the analysis code as run, with four mechanical edits, each verified by comparing
the abstract syntax tree against the working copy (logic identical):

1. Module docstrings and comments translated to English (all other code untouched).
2. `h3_panel_audit_unified.py`: figure output path changed from the internal Chinese build
   directory to `figures/` (created if absent).
3. `h1_sensitivity_battery.py`: one dead expression (`... if False else ""`) removed from a
   formatted string; the emitted value is unchanged.
4. `phylo_ed3.py`: two unused `import copy` statements removed.

`build_analysis_base_v1.py` restates the step that assembled the frozen analysis base from the
seasonised table; its column contract is stated in its docstring.

## Licence and citation

No licence has been chosen yet; add a `LICENSE` file before making the repository public. When the
paper is published, cite it instead of this repository. Contact: wchangx173@163.com.
