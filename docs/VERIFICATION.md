# Verification record — 25 November 2026

**Evidence:** console output and validated files supplied by the study author after execution on a Mac mini. This record is not an independent remote execution of R by the package reviewer.

- Environment reported in the local run: R **4.6.2** (Apple Silicon), macOS Tahoe 26.3; base R libraries only.
- Prepared input: 34 measurements × 793 absorbance features, with 25 pediatric participants, four new clean carrier strips, one archived strip, three normal-reference spectra and one gingivitis-reference spectrum.
- Author ran the five-stage analysis from the prepared input; fixed-carrier NMF errors at K = 0, 1, 2 and 3 reproduced **14.53645%, 6.08817%, 4.49977%, 3.25812%**.
- Both final figures were regenerated in the current repository's `R/results/` after correction of plotting-script paths.
- Terminal message: `PASS: 5/5 numeric checks, NMF summary present, 2 final PNG present`.
- Snapshot crosscheck: all 9 files under `outputs/` and both PNGs under `figures/` match the corresponding completed `R/results/` files in the author's uploaded archive byte-for-byte.
- The 5 numerical checks in `outputs/05_independent_numeric_checks.csv` are all `TRUE`; 11 public-input sanity checks in `checks/PUBLIC_INPUT_INDEPENDENT_PYTHON_QC.csv` are `TRUE`.
- The received ZIP included `R/results/` with restricted individual analytical files. These files were **excluded** from the distribution ZIP, along with macOS metadata. The distribution SHA-256 manifest was rebuilt after cleanup.
