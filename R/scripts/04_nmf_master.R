# Полностью новая NMF-ветка ОТ GCF_MASTER_34.csv (не от ранее подготовленного
# pediatric_X_active_absorbance.csv, происхождение которого не установлено).
# Отрицательные значения исходных спектров (~11%) не обнуляются: ко ВСЕМ
# образцам добавляется ОДИН И ТОТ ЖЕ положительный сдвиг -min(X) + epsilon.
# В факторизации фиксируются два известных спектра: чистый носитель и
# постоянный базисный вектор (учёт введённого аддитивного смещения).
# K = количество обучаемых спектральных компонент, K=0..3.
# Это математическая проверка состава вариации, а не выделение чистых веществ.

nmf_shift_bio <- max(0, -min(X_ped_bio)) + 1e-8
X_nmf_bio <- X_ped_bio + nmf_shift_bio
if (any(X_nmf_bio < 0)) stop("NMF: неверный неотрицательный сдвиг.")
n_bio <- nrow(X_nmf_bio)
p_bio <- ncol(X_nmf_bio)
norm_x_bio <- sqrt(sum(X_nmf_bio^2))
norm_strip_bio <- sqrt(sum(carrier_four_bio^2))
fixed_bio <- rbind(carrier_four_bio,
                   rep(norm_strip_bio / sqrt(p_bio), p_bio))
rownames(fixed_bio) <- c("fixed_carrier", "fixed_constant")

nmf_fit_bio <- function(K, seed, max_iter = 4000L) {
  if (K < 0L) stop("K должен быть неотрицательным.")
  set.seed(seed)
  if (K > 0L) {
    free <- matrix(runif(K * p_bio, min = 0.001, max = 0.05), nrow = K)
    H <- rbind(fixed_bio, free)
  } else H <- fixed_bio
  W <- matrix(runif(n_bio * (K + 2L), min = 0.01, max = 1), nrow = n_bio)
  if (K > 0L) for (j in 3L:(K + 2L)) {
    h_norm <- sqrt(sum(H[j, ]^2))
    W[, j] <- W[, j] * h_norm / norm_strip_bio
    H[j, ] <- H[j, ] * norm_strip_bio / h_norm
  }
  eps <- 1e-14
  old_error <- Inf
  small_change <- 0L
  done <- max_iter
  for (iteration in seq_len(max_iter)) {
    W <- W * pmax(X_nmf_bio %*% t(H), eps) /
      pmax(W %*% tcrossprod(H), eps)
    if (K > 0L) {
      free_rows <- 3L:(K + 2L)
      numerator <- t(W[, free_rows, drop = FALSE]) %*% X_nmf_bio
      denominator <- t(W[, free_rows, drop = FALSE]) %*% (W %*% H)
      H[free_rows, ] <- H[free_rows, , drop = FALSE] *
        pmax(numerator, eps) / pmax(denominator, eps)
      for (j in free_rows) {
        h_norm <- sqrt(sum(H[j, ]^2))
        if (!is.finite(h_norm) || h_norm <= eps) stop("NMF: вырожденная компонента.")
        W[, j] <- W[, j] * h_norm / norm_strip_bio
        H[j, ] <- H[j, ] * norm_strip_bio / h_norm
      }
    }
    if (iteration %% 100L == 0L) {
      current <- sqrt(sum((X_nmf_bio - W %*% H)^2)) / norm_x_bio
      if (!is.finite(current)) stop("NMF: нечисловая ошибка реконструкции.")
      if (abs(current - old_error) < 1e-9) small_change <- small_change + 1L
      else small_change <- 0L
      if (small_change >= 2L) { done <- iteration; break }
      old_error <- current
    }
  }
  list(W = W, H = H,
       error = sqrt(sum((X_nmf_bio - W %*% H)^2)) / norm_x_bio,
       iterations = done)
}

nmf_summary_rows_bio <- list()
nmf_fits_bio <- list()
nmf_seeds_bio <- 20261118L
for (K in 0:3) {
  starts <- if (K == 0L) 1L else 5L
  variants <- lapply(seq_len(starts), function(r)
    nmf_fit_bio(K, seed = nmf_seeds_bio + 1000L * K + r))
  errors <- vapply(variants, function(z) z$error, numeric(1))
  best <- variants[[which.min(errors)]]
  nmf_fits_bio[[as.character(K)]] <- best
  ch_idx <- which(wn_bio >= 2870 & wn_bio <= 3050)
  C <- X_nmf_bio[, ch_idx, drop = FALSE]
  P <- best$W %*% best$H[, ch_idx, drop = FALSE]
  max_cos <- if (K > 0L) {
    max(vapply(3L:(K + 2L), function(j)
      sum(best$H[j, ] * carrier_four_bio) /
        sqrt(sum(best$H[j, ]^2) * sum(carrier_four_bio^2)), numeric(1)))
  } else NA_real_
  nmf_summary_rows_bio[[K + 1L]] <- data.frame(
    K_free = K, fixed_components = 2L, starts = starts,
    error_full_percent = 100 * best$error,
    error_CH_percent = 100 * sqrt(sum((C - P)^2) / sum(C^2)),
    max_free_carrier_cosine = max_cos,
    best_iterations = best$iterations,
    worst_start_error_percent = 100 * max(errors),
    best_start_error_percent = 100 * min(errors))
  cat(sprintf("Новая NMF от master: K=%d, ошибка %.5f%%, стартов %d\n",
              K, 100 * best$error, starts))
}

nmf_summary_bio <- do.call(rbind, nmf_summary_rows_bio)
write.csv(nmf_summary_bio, file.path(result_dir, "04_NMF_master_summary.csv"),
          row.names = FALSE)
fit2 <- nmf_fits_bio[["2"]]
write.csv(data.frame(wavenumber = wn_bio, fixed_carrier = fit2$H[1L, ],
                     fixed_constant = fit2$H[2L, ],
                     C1 = fit2$H[3L, ], C2 = fit2$H[4L, ]),
          file.path(result_dir, "04_NMF_master_K2_components.csv"), row.names = FALSE)
nmf_coeff_bio <- data.frame(sample_id = meta_bio$sample_id,
                            coefficient_carrier = fit2$W[, 1L],
                            coefficient_constant = fit2$W[, 2L],
                            coefficient_C1 = fit2$W[, 3L],
                            coefficient_C2 = fit2$W[, 4L])
write.csv(nmf_coeff_bio,
          file.path(result_dir, "04_NMF_master_K2_coefficients_PRIVATE.csv"),
          row.names = FALSE)

# Ассоциация свободных коэффициентов с клиническим статусом и CH-PC1.
# Только описательный анализ после подгонки NMF без клинических меток.
nmf_association_bio <- do.call(rbind, lapply(c("coefficient_C1", "coefficient_C2"),
                                              function(cname) {
  z <- nmf_coeff_bio[[cname]]
  full <- exact_cliff_bio(z, case_bio)
  sel <- case_bio | age_control_bio
  age <- exact_cliff_bio(z[sel], case_bio[sel])
  data.frame(component = cname, cliff_all = full$effect,
             exact_p_all = full$exact_two_sided_p,
             cliff_age = age$effect, exact_p_age = age$exact_two_sided_p,
             spearman_CH_DPT4 = cor(z, pca_bio$EMSC_DPT4$CH$pc1_score,
                                     method = "spearman"))
}))
write.csv(nmf_association_bio, file.path(result_dir, "04_NMF_biology_exploratory.csv"),
          row.names = FALSE)
writeLines(c(
  paste0("Источник: ", basename(input_file)),
  sprintf("Доля исходных отрицательных спектральных значений: %.5f%%",
          100 * mean(X_ped_bio < 0)),
  sprintf("Один общий аддитивный сдвиг для всех спектров: %.12g", nmf_shift_bio),
  "Зафиксированы: средний спектр 4 чистых носителей и постоянный базис.",
  "Число обучаемых компонент: K=0,1,2,3.",
  "Алгоритм: евклидова NMF, мультипликативные обновления.",
  "Случайные зерна: 20261118 + 1000*K + номер запуска.",
  "Число запусков: K=0 один; K>0 пять. Не более 4000 итераций.",
  "Группы пародонтального статуса не использовались при подгонке NMF.",
  "Предыдущая NMF из отдельного active_absorbance.csv имеет иной вход и",
  "не должна численно сравниваться с новой NMF от GCF_MASTER_34.csv.",
  "",
  capture.output(sessionInfo())),
  file.path(result_dir, "04_NMF_method_and_R_versions.txt"))
