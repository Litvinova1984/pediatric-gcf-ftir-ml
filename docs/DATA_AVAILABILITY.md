# Data availability, ethics and scope

## Release candidate with a prepared public matrix

`data/public/GCF_MASTER_34_PUBLIC.csv` contains 34 prepared absorbance spectra × 793 wavenumbers, plus the minimal grouping variables needed for EMSC, unsupervised PCA, age-restricted exact permutation tests, fixed-carrier NMF and reconstruction of the final publication figures.

The 25 pediatric spectra stem from the same original participant measurements used in the authors' earlier Data study, **but this is a different analytical input**: that study released a 25 × 793 SNV-normalized pediatric matrix, whereas this study requires prepared absorbances for 25 children, four clean-strip controls, one archived strip, and four biological reference measurements. The two matrices must not be substituted for each other.

The public candidate removes source filenames, exact ages, sex, pubertal, endocrine and connective-tissue metadata. Pediatric records receive different analytical labels. It retains only the documented-periodontal-status grouping and the indicator of being in the 14–17-year comparison age band; both variables are indispensable for reproducing the paper's clinical comparisons. Figures and aggregates are also included.

## Ethics caveat for public publication

**Publication of the earlier Data spectra does not independently authorize linking previously published participant spectra to NEW periodontal-status information.** Even with reassigned IDs, spectral fingerprints may be linked across studies. Before public GitHub publication, the authorized research team must confirm that parental consent, ethics approval and local data governance cover disclosure of these specific participant-level spectral-and-clinical fields. If not, publish the clinical grouping only as aggregates and keep the complete input restricted; then state that full clinical reproducibility requires controlled access.

This ZIP is a *release candidate*. This note is not ethics approval.

## Computational scope

Run `Rscript run_all.R` from the repository root. Computation starts from the fixed prepared input matrix; upstream `.spc` and `.dpt` acquisition, preprocessing and its original conversion to this master table are not independently reproduced. The code uses base R. Computed individual scores and model coefficients go to `R/results/` and are Git-ignored. Both final PNG publication figures are created from computed tables.

## Licenses

MIT for R code; CC BY 4.0 for approved processed data, aggregate outputs, figures and documentation. This licensing statement does not authorize otherwise prohibited disclosure.
