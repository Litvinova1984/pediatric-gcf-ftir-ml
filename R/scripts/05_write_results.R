# Результаты в виде краткого проверяемого русскоязычного текста для рукописи.
# Текст генерируется из чисел ТЕКУЩЕГО запуска, а не вводится вручную.

take_bio <- function(df, carrier, region) {
  z <- df[df$carrier_reference == carrier & df$region == region, , drop = FALSE]
  if (nrow(z) != 1L) stop("Некорректный выбор модели в отчёте.")
  z[1L, ]
}
pc_ch <- take_bio(summary_bio, "EMSC_DPT4", "CH")
pc_protein <- take_bio(summary_bio, "EMSC_DPT4", "Protein")
main_ch <- take_bio(clinical_bio, "EMSC_DPT4", "CH")
main_protein <- take_bio(clinical_bio, "EMSC_DPT4", "Protein")
ch_stable <- take_bio(shape_qc_bio, "EMSC_DPT4", "CH")
cross_ch <- cross_shape_bio[cross_shape_bio$region == "CH", , drop = FALSE]
if (nrow(cross_ch) != 1L) stop("Нет межносительного CH-сопоставления.")
ext_g <- external_projection_bio[
  external_projection_bio$sample_group == "gingivitis_reference", , drop = FALSE]
ext_n <- external_projection_bio[
  external_projection_bio$sample_group == "normal_reference", , drop = FALSE]
if (nrow(ext_g) != 1L || nrow(ext_n) != 3L)
  stop("Некорректный набор внешних биологических референсов.")

paragraphs_bio <- c(
  "# Результаты биологической ветки: автоматически из текущего R-запуска",
  "",
  "**Статус:** вычислительная версия для сверки авторами; не окончательная редакция статьи.",
  "",
  "## Технический контроль и EMSC",
  "",
  sprintf(paste0("Исследованы 25 спектров детской десневой жидкости, ",
                 "четыре измерения чистого носителя и один архивный носитель. ",
                 "Минимальная корреляция между чистыми носителями составила r = %.6f; ",
                 "коэффициент вариации интегрального сигнала — %.3f%%."),
          min_carrier_cor, carrier_cv),
  "",
  "EMSC выполнена с основным шаблоном DPT4 и единым биологическим reference (средний спектр 25 детей). Как проверки чувствительности использованы проекции для DPT4 и архивного носителя с базисом из постоянной, линейной и квадратичной функций. Клинические метки в предобработке не использовались.",
  "",
  "## Спектральная вариация после EMSC",
  "",
  sprintf("В области 1540–1720 см⁻¹ PC1 объясняет %.2f%% дисперсии при DPT4; в области 2870–3050 см⁻¹ — %.2f%%.",
          pc_protein$PC1_variance_pct, pc_ch$PC1_variance_pct),
  "",
  "## Ассоциация с документированным пародонтальным статусом",
  "",
  sprintf(paste0("Для C–H-области по основной модели EMSC_DPT4 описательная ",
                 "величина эффекта Cliff's δ = %.3f, точное двустороннее ",
                 "перестановочное p = %.4f (3 против 22). ",
                 "При ограничении сравнительной группы шестью детьми 14–17 лет ",
                 "δ = %.3f и p = %.4f (3 против 6)."),
          main_ch$cliff_all, main_ch$exact_p_all,
          main_ch$cliff_age, main_ch$exact_p_age),
  "",
  sprintf(paste0("В белковой области 1540–1720 см⁻¹ соответствующая ",
                 "описательная величина эффекта составила δ = %.3f ",
                 "(точное p = %.4f); при возрастном сравнении δ = %.3f (p = %.4f)."),
          main_protein$cliff_all, main_protein$exact_p_all,
          main_protein$cliff_age, main_protein$exact_p_age),
  "",
  sprintf(paste0("Для C–H-области косинус направления группового ",
                 "спектрального контраста при полной и возрастной ",
                 "сравнительных группах составил %.3f. ",
                 "При исключении по одному случаю минимум составил %.3f. ",
                 "Контраст по EMSC_DPT4 был согласован с ",
                 "альтернативной проекционной коррекцией для архивного ",
                 "носителя (косинус %.3f)."),
          ch_stable$all_vs_age_contrast_cosine,
          ch_stable$leave_one_case_out_min_cosine,
          cross_ch$EMSC_vs_archive_projection_contrast_cosine),
  "",
  sprintf(paste0("При пассивной проекции четырёх независимых ",
                 "спектральных референсов в детское CH-PC1 пространство ",
                 "(EMSC с DPT4) образец G с документированным гингивитом ",
                 "получил координату %.3f SD; три контрольных образца N ",
                 "имели координаты в интервале от %.3f до %.3f SD. ",
                 "Это качественное дополнительное наблюдение, ",
                 "а не внешняя диагностическая валидация."),
          ext_g$CH_PC1_pediatric_SD, min(ext_n$CH_PC1_pediatric_SD),
          max(ext_n$CH_PC1_pediatric_SD)),
  "",
  sprintf(paste0("При этом корреляция индивидуальных координат CH-PC1 ",
                 "основной EMSC и архивной проекционной модели равна %.3f; ",
                 "устойчивость формы ",
                 "группового контраста не означает тождественность ",
                 "индивидуальных PCA-координат."),
          cross_ch$EMSC_vs_archive_projection_PC1_correlation),
  "",
  "## Биологическая интерпретация",
  "",
  "Результаты согласуются с наличием групповой вариации формы C–H-ассоциированных полос десневой жидкости, сопутствующей документированному состоянию пародонта. Эти полосы включают колебания C–H-связей, встречающихся в липидных и белковых структурах. Спектральный контраст сам по себе не идентифицирует конкретное вещество, его концентрацию или механизм воспаления.",
  "",
  "Из-за трёх случаев с документированной патологией, ассоциации пародонтального статуса с другими клиническими характеристиками и поискового характера сравнений нельзя заявлять диагностическую точность, специфичность или причинную связь. Приведённые p не скорректированы за множественные сравнения.",
  "")

if (exists("nmf_summary_bio", inherits = FALSE)) {
  e0 <- nmf_summary_bio[nmf_summary_bio$K_free == 0, ]
  e2 <- nmf_summary_bio[nmf_summary_bio$K_free == 2, ]
  paragraphs_bio <- c(paragraphs_bio,
    "## NMF как независимая проверка разложения спектра",
    "",
    sprintf(paste0("Для новой NMF, построенной непосредственно из ",
                   "GCF_MASTER_34.csv с единым положительным сдвигом, ",
                   "фиксированными спектром носителя и постоянным базисом, ",
                   "ошибка реконструкции равнялась %.3f%% при K=0 ",
                   "и %.3f%% при K=2."),
            e0$error_full_percent, e2$error_full_percent),
    "",
    "Эта NMF имеет иное входное преобразование, чем ранее сохранённая NMF из active_absorbance.csv; прямое сопоставление ошибок двух вариантов некорректно. Свободные компоненты не приравниваются к отдельным биохимическим веществам.",
    "")
}
paragraphs_bio <- c(paragraphs_bio,
  "## Воспроизводимость",
  "",
  "Все приведённые числа, групповые медианы, нагрузки и кривые построены из одной 34×793 матрицы при явных параметрах EMSC и фиксированных спектральных окнах. R-версия и вычислительный протокол сохраняются вместе с результатами.")
writeLines(paragraphs_bio, file.path(result_dir, "05_RESULTS_RU_AUTOGENERATED.md"),
           useBytes = TRUE)
# Контрольные числа получены НЕЗАВИСИМО другим вычислительным кодом из
# исходного CSV; не заменяют проверку запуска R на Mac.
independent_checks_bio <- data.frame(
  metric = c("EMSC DPT4 Protein PC1 %", "EMSC DPT4 CH PC1 %",
             "EMSC DPT4 CH two-sided permutation p",
             "CH full-vs-age contrast cosine",
             "CH EMSC-vs-archive-projection contrast cosine"),
  independent_reference = c(79.321984962, 76.534297475,
                            0.035652174, 0.972, 0.962020564),
  observed_R = c(pc_protein$PC1_variance_pct, pc_ch$PC1_variance_pct,
                 main_ch$exact_p_all,
                 ch_stable$all_vs_age_contrast_cosine,
                 cross_ch$EMSC_vs_archive_projection_contrast_cosine),
  absolute_tolerance = c(0.02, 0.02, 0.0006, 0.02, 0.02))
independent_checks_bio$pass <-
  abs(independent_checks_bio$independent_reference -
      independent_checks_bio$observed_R) <= independent_checks_bio$absolute_tolerance
write.csv(independent_checks_bio,
          file.path(result_dir, "05_independent_numeric_checks.csv"),
          row.names = FALSE)
if (!all(independent_checks_bio$pass))
  stop("Независимые контрольные величины не совпали: см. 05_independent_numeric_checks.csv")
writeLines(c("Сессия анализа GCF, biology_v2", capture.output(sessionInfo())),
           file.path(result_dir, "05_R_SESSION.txt"), useBytes = TRUE)
cat("Русский текст результатов сформирован из чисел текущего запуска.\n")
