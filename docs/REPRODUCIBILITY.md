# Reproducibility scope

## Starting point

The included file `data/public/GCF_MASTER_34_PUBLIC.csv` has 34 measurements and 793 prepared absorbance variables. It is distinct from the SNV-normalized pediatric matrix previously released with the authors' *Data* article.

From the root of the repository:

```bash
Rscript run_all.R
Rscript verify_local_results.R
```

The five stages are paper-strip/EMSC preprocessing, region-specific PCA, exploratory clinical association and stability checks, fixed-carrier NMF (K = 0–3), and results reporting. Two publication figures are generated afterwards. The project uses base R. Without macOS XQuartz/Cairo, optional PDF export is skipped, but high-resolution PNGs are produced.

`R/results/` is created automatically during execution and must remain excluded from publication: it contains individual scores, carrier coefficients and model matrices. `outputs/` and `figures/` hold the frozen aggregate output and two publication figures already checked by the research team. `checks/MAC_REFERENCE_2026-11-18.csv` provides numeric checkpoints, and `checks/PUBLIC_INPUT_INDEPENDENT_PYTHON_QC.csv` records 11 independent checks on the minimized public input.

## Validation on Mac mini

The author successfully ran the prepared-input numerical pipeline on an Apple Silicon Mac mini with R 4.6.2, then generated the final figures from `R/results/` within the same repository. After correcting the path selection and rerunning figure scripts, the verification command returned `PASS: 5/5 numeric checks, NMF summary present, 2 final PNG present`. The originally exported frozen aggregate files and figures match the corresponding files in the author's completed local `R/results/` byte-for-byte in the received archive.

The root `run_all.R` has been statically audited for repository-relative path construction, but an independent fresh-session execution of that **single command from this newly zipped archive** is not claimed. The user reports successful execution of the constituent stages and final verification.

## Boundary

The workflow reproduces numerical analyses and figures from the prepared absorbance matrix. It does **not** reproduce original instrument acquisition or upstream `.spc`/`.dpt` conversion to the 793-variable input. Privacy clearance for linking individual spectral rows to new clinical flags must be confirmed before the candidate is published.
