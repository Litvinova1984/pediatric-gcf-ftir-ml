# Сравнения клинических групп выполняются ПОСЛЕ EMSC/PCA.
# 22 без документированной патологии и 3 с документированной патологией;
# возрастная проверка: 6 без патологии (14–17 лет) против тех же 3 случаев.
# Результаты описательные: статус пародонта смешан с сопутствующими факторами.

cosine_bio <- function(a, b) {
  den <- sqrt(sum(a^2) * sum(b^2))
  if (den <= 0 || !is.finite(den)) return(NA_real_)
  sum(a * b) / den
}
cliff_bio <- function(a, b) mean(sign(outer(a, b, "-")))
exact_cliff_bio <- function(values, groups) {
  if (any(!is.finite(values)) || anyNA(groups)) stop("Недопустимые значения для permutation.")
  n <- length(values)
  n1 <- sum(groups)
  observed <- cliff_bio(values[groups], values[!groups])
  # Для 25 choose 3 = 2300 и 9 choose 3 = 84 возможен полный перебор.
  candidates <- combn(seq_len(n), n1)
  null <- apply(candidates, 2L, function(rows) {
    g <- rep(FALSE, n)
    g[rows] <- TRUE
    cliff_bio(values[g], values[!g])
  })
  list(effect = observed,
       exact_two_sided_p = mean(abs(null) >= abs(observed) - 1e-12),
       n_permutations = length(null))
}

shape_bio <- list()
comparison_rows_bio <- list()
shape_rows_bio <- list()
comparison_n <- 0L
shape_n <- 0L

for (nm in names(pca_bio)) {
  shape_bio[[nm]] <- list()
  for (region in names(regions_bio)) {
    z <- pca_bio[[nm]][[region]]
    ids_cases <- which(case_bio)
    ids_all_controls <- which(!case_bio)
    ids_age_controls <- which(age_control_bio)
    total <- exact_cliff_bio(z$pc1_score, case_bio)
    age_selection <- case_bio | age_control_bio
    age <- exact_cliff_bio(z$pc1_score[age_selection], case_bio[age_selection])
    loo <- vapply(ids_cases, function(i)
      cliff_bio(z$pc1_score[setdiff(ids_cases, i)],
                z$pc1_score[ids_all_controls]), numeric(1))

    comparison_n <- comparison_n + 1L
    comparison_rows_bio[[comparison_n]] <- data.frame(
      carrier_reference = nm, region = region,
      cases_n = length(ids_cases), controls_n = length(ids_all_controls),
      age_controls_n = length(ids_age_controls),
      median_case_PC1_SD = median(z$pc1_score[ids_cases]),
      median_control_PC1_SD = median(z$pc1_score[ids_all_controls]),
      median_age_control_PC1_SD = median(z$pc1_score[ids_age_controls]),
      cliff_all = total$effect, exact_p_all = total$exact_two_sided_p,
      cliff_age = age$effect, exact_p_age = age$exact_two_sided_p,
      loo_cliff_min = min(loo), loo_cliff_max = max(loo),
      permutations_all = total$n_permutations,
      permutations_age = age$n_permutations)

    M <- z$corrected
    m_case <- apply(M[ids_cases, , drop = FALSE], 2L, median)
    m_all <- apply(M[ids_all_controls, , drop = FALSE], 2L, median)
    m_age <- apply(M[ids_age_controls, , drop = FALSE], 2L, median)
    contrast_all <- m_case - m_all
    contrast_age <- m_case - m_age
    loo_cos <- vapply(ids_cases, function(i) {
      reduced <- apply(M[setdiff(ids_cases, i), , drop = FALSE], 2L, median) - m_all
      cosine_bio(contrast_all, reduced)
    }, numeric(1))
    shape_bio[[nm]][[region]] <- list(
      wn = z$wn, case = m_case, control_all = m_all, control_age = m_age,
      contrast_all = contrast_all, contrast_age = contrast_age,
      loo_direction_cos = loo_cos)
    shape_n <- shape_n + 1L
    shape_rows_bio[[shape_n]] <- data.frame(
      carrier_reference = nm, region = region,
      all_vs_age_contrast_cosine = cosine_bio(contrast_all, contrast_age),
      all_contrast_vs_PC1_loading_cosine = cosine_bio(contrast_all, z$pc1_loading),
      leave_one_case_out_min_cosine = min(loo_cos),
      leave_one_case_out_max_cosine = max(loo_cos),
      contrast_age_norm_over_all_norm =
        sqrt(sum(contrast_age^2) / sum(contrast_all^2)))

    write.csv(data.frame(wavenumber = z$wn, median_case = m_case,
                         median_control_all = m_all,
                         median_control_age = m_age,
                         contrast_all = contrast_all, contrast_age = contrast_age),
              file.path(result_dir, paste0("03_", nm, "_", region,
                                           "_group_profiles.csv")), row.names = FALSE)
  }
}

clinical_bio <- do.call(rbind, comparison_rows_bio)
shape_qc_bio <- do.call(rbind, shape_rows_bio)
cross_shape_bio <- do.call(rbind, lapply(names(regions_bio), function(region) {
  a <- shape_bio$EMSC_DPT4[[region]]
  p <- shape_bio$Projection_DPT4[[region]]
  b <- shape_bio$Projection_Archive[[region]]
  data.frame(region = region,
             EMSC_vs_archive_projection_contrast_cosine =
               cosine_bio(a$contrast_all, b$contrast_all),
             DPT4_vs_archive_projection_contrast_cosine =
               cosine_bio(p$contrast_all, b$contrast_all),
             EMSC_vs_archive_projection_age_cosine =
               cosine_bio(a$contrast_age, b$contrast_age),
             EMSC_vs_archive_projection_PC1_correlation =
               cor(pca_bio$EMSC_DPT4[[region]]$pc1_score,
                   pca_bio$Projection_Archive[[region]]$pc1_score))
}))

write.csv(clinical_bio, file.path(result_dir, "03_clinical_PC1_exact_permutations.csv"),
          row.names = FALSE)
write.csv(shape_qc_bio, file.path(result_dir, "03_shape_stability.csv"),
          row.names = FALSE)
write.csv(cross_shape_bio, file.path(result_dir, "03_cross_carrier_sensitivity.csv"),
          row.names = FALSE)

png(file.path(result_dir, "Figure_CH_biological_contrast.png"),
    width = 1900, height = 920, res = 180)
par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3, 1))
for (nm in c("EMSC_DPT4", "Projection_Archive")) {
  z <- shape_bio[[nm]]$CH
  ylim <- range(c(z$contrast_all, z$contrast_age)) * 1.1
  plot(z$wn, z$contrast_all, type = "l", lwd = 2.5, col = "#345D82",
       xlim = rev(range(z$wn)), ylim = ylim,
       main = paste("CH contrast, carrier:", nm),
       xlab = "Wavenumber (cm-1)", ylab = "Median spectral contrast")
  lines(z$wn, z$contrast_age, col = "#C25B41", lwd = 2, lty = 2)
  abline(h = 0, lty = 3, col = "grey60")
  legend("topright", legend = c("3 vs 22", "3 vs age-matched 6"),
         col = c("#345D82", "#C25B41"), lwd = 2, lty = c(1, 2), bty = "n", cex = 0.78)
}
dev.off()

cat("Сравнение клинических групп (точные перестановки):\n")
print(clinical_bio[, c("carrier_reference", "region", "cliff_all", "exact_p_all",
                       "cliff_age", "exact_p_age")], row.names = FALSE, digits = 5)
cat("Форма межгруппового спектрального контраста и устойчивость:\n")
print(shape_qc_bio, row.names = FALSE, digits = 4)
print(cross_shape_bio, row.names = FALSE, digits = 4)
