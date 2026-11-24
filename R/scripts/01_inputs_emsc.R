# Входные данные, технический контроль и EMSC, заданная аналитически.
# Спектры остаются неизменными; EMSC оценивает форму после учёта фона,
# мультипликативной вариации и известного спектра чистого носителя.

if (!file.exists(input_file)) stop("Нет подготовленной спектральной матрицы: ", input_file)
master_bio <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
features <- grep("^wn_[0-9]+(\\.[0-9]+)?$", names(master_bio), value = TRUE)
needed <- c("sample_id", "sample_group", "periodontal_status")
if ("age_14_17" %in% names(master_bio)) {
  needed <- c(needed, "age_14_17")
} else if ("age" %in% names(master_bio)) {
  needed <- c(needed, "age")
} else {
  stop("Не найден возрастной признак 14–17 лет (или возраст во входе с ограниченным доступом).")
}
if (nrow(master_bio) != 34L || length(features) != 793L ||
    !all(needed %in% names(master_bio))) stop("Неверная структура исходной матрицы.")

wn_bio <- as.numeric(sub("^wn_", "", features))
if (anyNA(wn_bio) || any(diff(wn_bio) >= 0) || anyDuplicated(wn_bio))
  stop("Волновые числа должны идти по убыванию без повторов.")
X_full_bio <- as.matrix(master_bio[, features, drop = FALSE])
storage.mode(X_full_bio) <- "double"
if (any(!is.finite(X_full_bio))) stop("В спектрах есть NA/NaN/Inf.")

ped_positions <- which(master_bio$sample_group == "pediatric")
# Сохраняем зафиксированный порядок строк публичного входа; никаких числовых ID пациентов.
clean_positions <- which(master_bio$sample_group == "clean_strip")
archive_positions <- which(master_bio$sample_group == "archive_strip")
if (length(ped_positions) != 25L || length(clean_positions) != 4L ||
    length(archive_positions) != 1L ||
    anyDuplicated(as.character(master_bio$sample_id[ped_positions])))
  stop("Неверный состав группы детей или носителей.")
X_ped_bio <- X_full_bio[ped_positions, , drop = FALSE]
meta_bio <- master_bio[ped_positions, needed, drop = FALSE]
case_bio <- meta_bio$periodontal_status == "Documented pathology"
# Проверяем только размер клинической группы; не раскрываем сопоставление
# аналитических идентификаторов с диагнозом в публичном исходном коде.
if (anyNA(case_bio) || sum(case_bio) != 3L)
  stop("Ожидаются три наблюдения с документированной патологией пародонта.")
if ("age_14_17" %in% names(meta_bio)) {
  eligible_age_bio <- as.logical(meta_bio$age_14_17)
} else {
  eligible_age_bio <- !is.na(meta_bio$age) & meta_bio$age >= 14 & meta_bio$age <= 17
}
if (anyNA(eligible_age_bio)) stop("Пропущена принадлежность к возрастной группе 14–17 лет.")
age_control_bio <- !case_bio & eligible_age_bio
if (sum(age_control_bio) != 6L || !all(eligible_age_bio[case_bio]))
  stop("Возрастной контроль должен содержать 6 детей 14–17 лет; все 3 случая также 14–17 лет.")

carrier_four_bio <- colMeans(X_full_bio[clean_positions, , drop = FALSE])
carrier_archive_bio <- drop(X_full_bio[archive_positions, ])
carrier_list_bio <- list(DPT4 = carrier_four_bio, Archive = carrier_archive_bio)
carrier_cor_bio <- cor(t(X_full_bio[clean_positions, , drop = FALSE]))
min_carrier_cor <- min(carrier_cor_bio[upper.tri(carrier_cor_bio)])
carrier_integrals <- rowSums(X_full_bio[clean_positions, , drop = FALSE])
carrier_cv <- 100 * sd(carrier_integrals) / mean(carrier_integrals)
if (abs(min_carrier_cor - 0.999988991547) > 1e-6 ||
    abs(carrier_cv - 0.818884714) > 0.0005)
  stop("Чистые носители не совпадают с независимыми контрольными значениями.")
write.csv(data.frame(n_children = 25L, n_documented_periodontal = sum(case_bio),
                     n_age_control = sum(age_control_bio),
                     n_clean_carriers = 4L,
                     min_clean_carrier_correlation = min_carrier_cor,
                     cv_integrated_carrier_percent = carrier_cv),
          file.path(result_dir, "01_carrier_qc.csv"), row.names = FALSE)

# EMSC с ОСНОВНЫМ шаблоном DPT4: x_i(v) = a_i + l_i*t(v) + q_i*t(v)^2 + b_i*reference(v)
#            + c_i*carrier(v) + epsilon_i(v).
# corrected_i(v) = [x_i(v) - a_i - l_i*t(v) - q_i*t(v)^2
#                   - c_i*carrier(v)] / b_i.
# Reference: средний спектр 25 детей БЕЗ использования клинических меток.
# Архивный шаблон не используется в нормировке на b_i: при подстановке
# archive в полную EMSC коэффициент b_i у части детей отрицателен.
# Вместо такой невалидной нормировки проверяем устойчивость через
# определённую ниже проекцию на nuisance-базис (1,t,t²,carrier).
t_bio <- (wn_bio - mean(wn_bio)) / diff(range(wn_bio))
reference_bio <- colMeans(X_ped_bio)
fit_emsc_bio <- function(M, S) {
  design <- cbind(intercept = rep(1, length(wn_bio)),
                  linear = t_bio, quadratic = t_bio^2,
                  pediatric_reference = reference_bio, carrier = S)
  q <- qr(design, tol = 1e-10)
  if (q$rank != ncol(design)) stop("EMSC: зависимые столбцы дизайна.")
  beta <- t(qr.coef(q, t(M)))
  colnames(beta) <- colnames(design)
  if (any(!is.finite(beta)) || any(beta[, "pediatric_reference"] <= 0.1))
    stop("EMSC: неустойчивый или неположительный множитель reference.")
  subtracted <- outer(beta[, "intercept"], design[, "intercept"]) +
    outer(beta[, "linear"], design[, "linear"]) +
    outer(beta[, "quadratic"], design[, "quadratic"]) +
    outer(beta[, "carrier"], design[, "carrier"])
  corrected <- sweep(M - subtracted, 1L, beta[, "pediatric_reference"], "/")
  if (any(!is.finite(corrected))) stop("EMSC: недопустимые скорректированные значения.")
  list(corrected = corrected, coefficients = beta,
       condition_number = kappa(design))
}

project_nuisance_bio <- function(M, S) {
  design <- cbind(intercept = rep(1, length(wn_bio)),
                  linear = t_bio, quadratic = t_bio^2, carrier = S)
  q <- qr(design, tol = 1e-10)
  if (q$rank != ncol(design)) stop("Проекция: зависимые столбцы дизайна.")
  coeff <- t(qr.coef(q, t(M)))
  colnames(coeff) <- colnames(design)
  corrected <- M - coeff %*% t(design)
  if (any(!is.finite(corrected))) stop("Проекция: неконечный результат.")
  list(corrected = corrected, coefficients = coeff,
       condition_number = kappa(design))
}

preprocessed_bio <- list(
  EMSC_DPT4 = fit_emsc_bio(X_ped_bio, carrier_four_bio),
  Projection_DPT4 = project_nuisance_bio(X_ped_bio, carrier_four_bio),
  Projection_Archive = project_nuisance_bio(X_ped_bio, carrier_archive_bio))
for (nm in names(preprocessed_bio)) {
  write.csv(cbind(sample_id = meta_bio$sample_id,
                  as.data.frame(preprocessed_bio[[nm]]$coefficients)),
            file.path(result_dir, paste0("01_", nm, "_coefficients_PRIVATE.csv")),
            row.names = FALSE)
}
saveRDS(preprocessed_bio, file.path(result_dir, "01_corrected_matrices_PRIVATE.rds"))
cat(sprintf("4 носителя: min r=%.9f, CV=%.4f%%. Дети: 22+3; возрастной контроль: 6.\n",
            min_carrier_cor, carrier_cv))
cat(sprintf(paste0("EMSC DPT4: число обусловленности %.2f; ",
                   "проекция DPT4 %.2f, архивная проекция %.2f.\n"),
            preprocessed_bio$EMSC_DPT4$condition_number,
            preprocessed_bio$Projection_DPT4$condition_number,
            preprocessed_bio$Projection_Archive$condition_number))
