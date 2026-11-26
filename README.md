# Corrected AvR diagnostic component

Owner: Changxian Wu (CW), College of Agriculture and Animal Husbandry Science and Technology, Ningxia Vocational and Technical College for Minorities, Wuzhong 751100, Ningxia, China. Contact: wchangx173@163.com.

This two-file code package reproduces the corrected diagnostic oligonucleotide component from supplied, previously extracted RSV gene windows and ENA metadata. It is **component reproduction, not reproduction of the entire study**. It does not download sequences, reconstruct the input windows, assign subtypes, fit exposure models, infer phylogenies, or regenerate submission documents or figures. No sequence or metadata datasets or generated analysis tables are bundled. This local package does not establish that code or inputs have been publicly uploaded or deposited.

## Requirements

Python 3.10 or later, NumPy, and the third-party `regex` package. No pandas, Office software, aligner, or internet connection is required to execute the analysis once the dependencies and inputs are available.

```bash
python -m pip install numpy regex
```

Dependencies are not version-pinned in this two-file package. Each analysis writes its actual Python, NumPy, and regex versions, script hash, input hashes, and output hashes to `provenance.json`.

## Exact inputs

Supply these four files directly inside the directory passed to `--input-dir`, preserving the filenames:

| Filename | Required contents and selection |
|---|---|
| `ena_rsv_meta_matrix_sp11250_2026-10-27.tsv` | Tab-separated ENA metadata with header columns `accession`, `collection_date`, and `first_public`; additional columns are allowed. Accessions must be unique after removing version suffixes. |
| `win_M_full.fasta` | Existing orientation-normalized M-gene windows; script selects explicit subtype A headers. |
| `win_L_B.fasta` | Existing orientation-normalized, **subtype-B-only** L-gene windows. Untagged headers are accepted on this assumption; an explicit non-B tag is excluded. |
| `win_n_full.fasta` | Existing orientation-normalized N-gene windows; script selects explicit subtype B headers. Filename case is significant on case-sensitive filesystems. |

M and N headers must start with an accession and identify subtype, for example `>ACCESSION.1 A`, `>ACCESSION.1 B`, or `>ACCESSION.1 subtype=B`. For L, `>ACCESSION.1` is accepted only because the caller supplies an already subtype-B-only file. Sequences may span lines. They are uppercased; non-ACGT symbols remain ambiguous, not wildcards. A duplicate accession within a selected gene cohort is rejected. Metadata is joined by accession without its version suffix; records without matching metadata have unresolved dates and are excluded from annual analyses, with counts retained in provenance outputs.

These are prepared intermediate inputs. Obtaining a current metadata download or independently extracting new windows does not guarantee the same cohort or results as a frozen input set. Archive the exact inputs separately with their provenance if numerical reproduction is required.

## Running

Run the script from a trusted working directory using absolute input/output paths. Keep input data separate from code.

```bash
python -I /absolute/path/code_release_avr_revision/audit_diagnostics.py --self-test-only
python -I /absolute/path/code_release_avr_revision/audit_diagnostics.py --input-dir /absolute/path/inputs --output-dir /absolute/path/new_results
```

`--selftest` is an alias for `--self-test-only`. Self-tests require no input data, create no result directory, and print the number passed. Analysis runs self-tests first. The module has a main guard: importing it does not start analysis or write results.

`--input-dir` and `--output-dir` are required for analysis. The result directory must be new or empty, outside the release directory, and separate from and not nested with the input directory. Existing results are not overwritten. The CLI reads and hashes only the four supplied inputs; it has no dependency on legacy analysis files, submission documents, or submission figures. Interrupted runs can leave partial results; use another empty directory to rerun. Input hashes are checked again after a successful run.

## Dates and denominators

Analysis covers calendar years 2015–2026; the detailed target year is 2026. These constants are fixed in the script. Collection dates are parsed first with calendar validation: `yyyy-mm-dd`, `yyyy-mm`, `yyyy`, `yyyy-Mon`, `Mon-yyyy`, and `dd-Mon-yyyy`, with English month abbreviations. Same-year ranges resolve a year; cross-year ranges do not specify one collection year. Missing or unparseable collection dates use a parsed `first_public` year only as an explicitly labelled fallback. **Publication year is not collection year.** Records with neither resolvable date are excluded from annual cohorts. All date sources, formats, fallback reasons, and counts are exported. A separate collection-date-only sensitivity excludes every publication fallback.

The denominator `n` is the number of distinct accession records with a supplied gene/subtype window and the relevant year. It is fixed for every arm in that cohort/year. Ambiguous, too-short, or unmatched windows remain in `n`; records without that gene window are outside that cohort. Thus `n` is not the number of all available public genomes. Subtypes, window extraction, and orientation normalization are accepted as prepared-input assumptions and are not independently inferred by this component.

## Matching and positions

The 17 oligonucleotide arms are embedded in `ARMS`. The script searches both orientations for all overlapping equal-length, entirely A/C/G/T substrings within two substitutions, using `regex` constraints `{s<=2,i<=0,d<=0}` and independent Hamming checks. It retains all sites with the globally minimum distance. Insertions and deletions are not allowed. Ambiguous substrings cannot qualify as exact or approximate matches; an otherwise ambiguous window can still contain a qualifying canonical site.

`p_exact`, `p_le1`, and `p_le2` are cumulative per-record fractions on the fixed denominator. No canonical hit within two substitutions can reflect ambiguity, missing sites, indels, or sequence divergence; it is not proof of three or more biological substitutions. Exact distances above two are not estimated.

Positions are one-based from the original oligo's 3′ end: for zero-based oriented-query index `i`, forward position is `L-i`, reverse position is `i+1`. Output includes all `frac_pos_-1` through `frac_pos_-5`. `p_3prox` records a mismatch anywhere in those five bases. Internal positions exclude both terminal five-base regions; `frac_pos_internal` and `frac_pos_5p` count a record once per category, not once per mismatch. Categories are marginal and need not sum to an overall mismatch fraction.

Primary position fractions use the union across tied best hits; `*_all_best_hits` requires the category in every tied hit. Positional disagreements are counted, and the per-record audit retains all hits plus a deterministic representative. Positions are observable only among canonical best hits within two substitutions. Unmatched records contribute no observed position; these frequencies are not complete positional mismatch prevalence. Whole-window best-site search does not establish a unique biological binding site. A composite primer arm's functional role must qualify positional interpretation.

Annual exact-compatibility drift is descriptive weighted least squares of percentage versus **calendar year** during 2018–2026, weighted by annual `n`. Included years require `n >= 30`, at least three years spanning at least three years. Both collection-first and collection-only slopes are exported. Publication fallbacks and changes in sequence sampling can affect slopes; no causal or assay-performance conclusion follows.

## Generated outputs

- `per_arm_2026.tsv`: detailed collection-first results, including fallback counts, every positional marginal, ambiguity/unmatched counts, and tied-hit diagnostics. Fractions are stored on a 0–1 scale.
- `per_arm_2026_collection_only.tsv`: target-year sensitivity without publication fallback.
- `per_arm_by_year.tsv`: annual compatibility and position results for both date bases.
- `exact_drift_2018_2026.tsv`: descriptive annual slopes and years/records used.
- `date_source_counts.tsv` and `cohort_record_dates.tsv`: date provenance and cohort exclusions.
- `per_record_audit.tsv.gz`: per-record/per-arm match evidence, source coordinates, and all best hits.
- `synthetic_tests.json`, `methods_and_limitations.json`, and `provenance.json`: validation, definitions, and reproducibility metadata.

The built-in tests cover exact and one/two/three substitutions, both orientations, terminal indexing, internal/5′ record marginals, indel and ambiguity handling, date parsing/fallbacks, tied sites, and comparison with exhaustive Hamming scans. These checks validate the component logic; they do not validate upstream extraction or demonstrate wet-lab performance. No wet-work results are generated or implied.
