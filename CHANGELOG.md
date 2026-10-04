# Changelog

## v1.0.1 (2026-10-04)

### Fixed
- **Crash on scientific-notation coordinates.** TRASH writes its CSV tables with R `write.csv`,
  which prints round coordinates in scientific notation (e.g. `2.8e+07`). pandas then read the
  whole `start`/`end` column as float, `candidate_arrays.bed` was written with `.0` suffixes, and
  monomer extraction (`modules/monomer_core.py`, both the standard and `--refined` modes) failed
  with `ValueError` on `int()`. Coordinates are now cast to integers when the TRASH tables are
  read (`modules/array_parser.py`, `modules/monomer_core.py`). Whether v1.0.0 crashed depended on
  the data (a single such coordinate is enough); when it did not crash, the only effect was a `.0`
  in monomer FASTA headers — extracted sequences are identical.

### Changed
- **TRASH output table selection (behavior change).** v1.0.0 took the *largest* file matching
  `*arrays.csv`, which in practice was TRASH's temporary `*_classarrays.csv` (an intermediate
  table that still contains arrays with zero mapped repeats). v1.0.1 reads TRASH's final
  `<genome>_arrays.csv` by name (excluding `*_no_repeats_arrays.csv`), requires exactly one final
  arrays table and one `*repeats_with_seq.csv` under the TRASH directory, and checks that the
  arrays table has the `repeats_number` column (i.e. TRASH finished). Otherwise AutoCen stops with
  an error listing what it found.
  - Consequence: candidate arrays — and therefore satellite families — can differ from v1.0.0
    results on the same TRASH output. Results produced with v1.0.0 remain reproducible with
    v1.0.0 (version DOI `10.5281/zenodo.21475673`).
  - A `--trash_dir` directory must now hold exactly one finished TRASH run.

## v1.0.0 (2026-07-21)

Initial release.
