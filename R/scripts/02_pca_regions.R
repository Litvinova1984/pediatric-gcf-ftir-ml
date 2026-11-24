# Без клинических меток: PCA скорректированных спектров в двух фиксированных
# областях. Варианты DPT4/Archive обучаются независимо; знак PC1 определяется
# только по нагрузкам, а не по различию клинических групп. EMSC_DPT4 —
# основная модель; Projection_DPT4/Archive — проверки чувствительности.

regions_bio <- list(Protein = c(1540, 1720), CH = c(2870, 3050))
fit_region_bio <- function(M, bounds) {
  idx <- which(wn_bio >= bounds[1L] & wn_bio <= bounds[2L])
  W <- M[, idx, drop = FALSE]
  if (length(idx) < 10L || nrow(W) < 4L) stop("Неожиданный размер спектрального окна.")
  p <- prcomp(W, center = TRUE, scale. = FALSE)
  loading <- p$rotation[, 1L]
  orientation <- if (loading[which.max(abs(loading))] < 0) -1 else 1
  raw_score <- orientation * p$x[, 1L]
  score_sd <- sd(raw_score)
  if (!is.finite(score_sd) || score_sd <= 0) stop("PC1 с нулевой дисперсией.")
  list(indices = idx, wn = wn_bio[idx], corrected = W,
       pc1_score = as.numeric(raw_score / score_sd),
       pc1_loading = orientation * loading,
       pc1_explained_pct = 100 * p$sdev[1L]^2 / sum(p$sdev^2),
       raw_score_sd = score_sd,
       pediatric_spectral_mean = colMeans(W))
}

pca_bio <- lapply(names(preprocessed_bio), function(nm) {
  lapply(regions_bio, function(bounds)
    fit_region_bio(preprocessed_bio[[nm]]$corrected, bounds))
})
names(pca_bio) <- names(preprocessed_bio)

summary_bio <- do.call(rbind, lapply(names(pca_bio), function(nm) {
  do.call(rbind, lapply(names(regions_bio), function(region) {
    x <- pca_bio[[nm]][[region]]
    data.frame(carrier_reference = nm, region = region,
               lower_cm1 = min(x$wn), upper_cm1 = max(x$wn),
               n_wavenumbers = length(x$wn),
               PC1_variance_pct = x$pc1_explained_pct,
               PC1_sd_before_standardization = x$raw_score_sd)
  }))
}))
write.csv(summary_bio, file.path(result_dir, "02_pca_region_summary.csv"),
          row.names = FALSE)

scores_bio <- data.frame(sample_id = meta_bio$sample_id,
                         periodontal_status = meta_bio$periodontal_status,
                         age_14_17 = if ("age_14_17" %in% names(meta_bio))
                           meta_bio$age_14_17 else (meta_bio$age >= 14 & meta_bio$age <= 17),
                         stringsAsFactors = FALSE)
for (nm in names(pca_bio)) for (region in names(regions_bio)) {
  scores_bio[[paste(nm, region, "PC1", sep = "_")]] <-
    pca_bio[[nm]][[region]]$pc1_score
}
write.csv(scores_bio, file.path(result_dir, "02_pca_scores_PRIVATE.csv"),
          row.names = FALSE)

# Четыре НЕЗАВИСИМЫХ биологических референса (N, N-1, N-2, G)
# проецируются ПОСЛЕ фиксации EMSC reference и детской PCA, без refit.
# Для внешней проверки используем основной шаблон DPT4, где EMSC beta>0.
external_positions <- which(master_bio$sample_group %in%
                              c("normal_reference", "gingivitis_reference"))
if (length(external_positions) != 4L) stop("Ожидаются четыре внешних спектра.")
external_matrix <- X_full_bio[external_positions, , drop = FALSE]
external_fit <- fit_emsc_bio(external_matrix, carrier_four_bio)
external_projection_bio <- data.frame(
  sample_id = master_bio$sample_id[external_positions],
  sample_group = master_bio$sample_group[external_positions],
  EMSC_reference_beta = external_fit$coefficients[, "pediatric_reference"],
  stringsAsFactors = FALSE)
for (region in names(regions_bio)) {
  z <- pca_bio$EMSC_DPT4[[region]]
  ext <- external_fit$corrected[, z$indices, drop = FALSE]
  centered <- sweep(ext, 2L, z$pediatric_spectral_mean, "-")
  external_projection_bio[[paste0(region, "_PC1_pediatric_SD")]] <-
    as.numeric(centered %*% z$pc1_loading / z$raw_score_sd)
}
write.csv(external_projection_bio,
          file.path(result_dir, "02_external_reference_projection.csv"),
          row.names = FALSE)
saveRDS(pca_bio, file.path(result_dir, "02_pca_models_PRIVATE.rds"))
cat("PCA рассчитана с нуля в двух областях для обоих носителей.\n")
print(summary_bio[, c("carrier_reference", "region", "PC1_variance_pct")],
      row.names = FALSE, digits = 5)
