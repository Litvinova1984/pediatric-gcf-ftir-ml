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



## Boundary

The workflow reproduces numerical analyses and figures from the prepared absorbance matrix. It does **not** reproduce original instrument acquisition or upstream `.spc`/`.dpt` conversion to the 793-variable input.  
