# Per-arm sequence compatibility of published RSV diagnostic oligonucleotides

Analysis code for the manuscript:

> **Per-arm sequence compatibility of 17 published respiratory syncytial virus diagnostic
> oligonucleotides with contemporary genomes: a reproducible in silico audit**

Submitted to the *Journal of Virological Methods*. Sole author and corresponding author:
Changxian Wu, College of Agriculture and Animal Husbandry Science and Technology, Ningxia
Vocational and Technical College for Minorities, Wuzhong 751100, Ningxia, China
(wchangx173@163.com; ORCID 0000-0003-2486-3380).

This repository contains analysis code only: no sequence data, no derived tables, no credentials.
Every sequence analysed is public and is re-downloaded by the scripts in `scripts/`.

---

## Run the audit component (the part the paper reports)

`audit_diagnostics.py` (repository root) reproduces the per-arm audit from four prepared inputs and
writes a hashed provenance record. It is a single file with no dependency on the rest of the
repository.

```bash
python -m pip install numpy regex
python -I audit_diagnostics.py --self-test-only        # 130 synthetic checks, no data needed
python -I audit_diagnostics.py --input-dir /absolute/path/inputs \
                              --output-dir /absolute/path/new_results
```

Inputs, all placed directly in `--input-dir`:

| File | Required content |
|---|---|
| `ena_rsv_meta_matrix_sp11250_2026-10-27.tsv` | ENA metadata with `accession`, `collection_date`, `first_public` |
| `win_M_full.fasta` | orientation-normalised M-gene windows, subtype A headers |
| `win_L_B.fasta` | orientation-normalised L-gene windows, subtype B only |
| `win_n_full.fasta` | orientation-normalised N-gene windows, subtype B headers |

Outputs are the per-arm 2026 table, the annual tables, the drift table, the date-provenance
counts, the per-record audit, and `provenance.json` (script hash, input hashes, output hashes,
interpreter and library versions). Subtypes, window extraction and orientation normalisation are
accepted as prepared-input assumptions; the script does not infer them.

Method in one sentence: each of the 17 published oligonucleotide arms is searched, in both
orientations, against every window of its own gene and subtype, allowing substitutions only, and
the mismatches of the best hit are resolved as positions measured from the arm's 3′ end.

---

## The rest of the pipeline

| Layer | Question | Script |
|---|---|---|
| Retrieval | all public human RSV nucleotide records | `download_ed3_nuccore_fasta.py` |
| Windows | F-region windows by 12-mer seed density | `extract_f_windows.py` |
| Windows | N-region windows | `extract_n_windows.py` |
| Alignment input | reference frame and orientation normalisation | `prep_mafft_inputs.py` |
| Variant calling | codon-aware substitutions in the reference frame | `call_subs_from_aln.py` |
| Trajectories | substitution frequency by year and country group | `traj_site0_by_year.py` |
| Seasons | collection season and the analysis base | `seasonize_and_h1base.py`, `build_analysis_base_v1.py` |
| Marker share | B.D.E.1 marker-set share by collection season | `make_season_lineage_table.py` |
| Audit (per arm) | reference-free exact and substitution-only search | `h3_panel_audit_unified.py` |
| Audit (position) | position-resolved mismatch profiles | `h3_panel_position_profile.py` |
| Audit (counts) | exact / one- / two-substitution distributions | `h3_panel_mismatch_dist.py` |
| Audit (pilot) | first-pass probe audit and sanity checks | `h3_probe_audit_pilot.py` |
| Phylogeny | F+G tree input and tree running | `phylo_ed3.py` |
| Context model | exploratory programme-exposure model (not a result of this paper) | `h1_frozen_outcome_run.py`, `h1_sensitivity_battery.py` |

Requirements are listed in `requirements.txt`. Reference sequences used for the F-region frame are
in `scripts/refseq_rsv_fasta_2026-10-28.txt`.

---

## Scope and limits

- The repository reproduces the audit; it does not reproduce the deposit-level filtering that
  produced the frozen window files, which are prepared intermediates described in the paper.
- No sequence data or derived tables are deposited here. Obtaining a current download does not
  guarantee the same cohort as the frozen one used in the paper.
- The exploratory programme-exposure analysis in the supplement is context only; it is
  inconclusive and the scripts are included for completeness, not as a reported result.
