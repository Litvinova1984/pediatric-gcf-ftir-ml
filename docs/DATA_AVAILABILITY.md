# Data availability, ethics and scope

## Release candidate with a prepared public matrix

`data/public/GCF_MASTER_34_PUBLIC.csv` contains 34 prepared absorbance spectra × 793 wavenumbers, plus the minimal grouping variables needed for EMSC, unsupervised PCA, age-restricted exact permutation tests, fixed-carrier NMF and reconstruction of the final publication figures.

The public candidate removes source filenames, exact ages, sex, pubertal, endocrine and connective-tissue metadata. Pediatric records receive different analytical labels. It retains only the documented-periodontal-status grouping and the indicator of being in the 14–17-year comparison age band; both variables are indispensable for reproducing the paper's clinical comparisons. Figures and aggregates are also included.

## Computational scope

Run `Rscript run_all.R` from the repository root. Computation starts from the fixed prepared input matrix; upstream `.spc` and `.dpt` acquisition, preprocessing and its original conversion to this master table are not independently reproduced. The code uses base R. Computed individual scores and model coefficients go to `R/results/` and are Git-ignored. Both final PNG publication figures are created from computed tables.

## Licenses

MIT for R code; CC BY 4.0 for approved processed data, aggregate outputs, figures and documentation. This licensing statement does not authorize otherwise prohibited disclosure.
