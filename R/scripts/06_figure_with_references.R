# Рисунок 2. Русская публикационная версия с разметкой по правкам Павла.
# Изображение строится из сохранённых результатов EMSC/PCA, НЕ из готового PNG.
# Статистика, порядок пациентов, коэффициенты EMSC и оси PCA не пересчитываются.
# Декоративная разметка является частью кода: скобка ν_as(CH2)
# и подсветка 2915–2935 см^-1 (точные границы зафиксированы ниже).
#
# Запуск: source(file.choose()), выбрав этот скрипт из biology_v2/scripts/.
# Входы должны находиться в biology_v2/results/.
# Только base R; печать PNG и, если доступен Cairo, PDF с кириллицей.

figure_bio_dir <- function() {
  explicit_root <- getOption("gcf.figure.biology_dir", NULL)
  if (!is.null(explicit_root))
    return(normalizePath(explicit_root, mustWork = TRUE))
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  paths <- paths[!is.na(paths) & basename(paths) == "06_figure_with_references.R"]
  if (length(paths)) {
    script <- normalizePath(tail(paths, 1L), mustWork = TRUE)
    return(normalizePath(file.path(dirname(script), ".."), mustWork = TRUE))
  }
  # Never read an old global biology_dir left by a previous RStudio run.
  stop("Запустите 06_figure_with_references.R через source(file.choose()).")
}

make_bio_figure_with_refs <- function() {
  root <- figure_bio_dir()
  results_dir <- file.path(root, "results")
  input_files <- file.path(results_dir, c(
    "03_EMSC_DPT4_CH_group_profiles.csv",
    "02_pca_scores_PRIVATE.csv",
    "02_external_reference_projection.csv"
  ))
  if (any(!file.exists(input_files)))
    stop("Нет нужных CSV в biology_v2/results/. Сначала выполните 00_run_biology.R.\n",
         paste(input_files[!file.exists(input_files)], collapse = "\n"))

  profiles <- read.csv(input_files[1], check.names = FALSE)
  scores <- read.csv(input_files[2], check.names = FALSE)
  external <- read.csv(input_files[3], check.names = FALSE)

  if (!all(c("wavenumber", "contrast_all", "contrast_age") %in% names(profiles)) ||
      !all(c("periodontal_status", "EMSC_DPT4_CH_PC1") %in% names(scores)) ||
      !all(c("sample_id", "sample_group", "CH_PC1_pediatric_SD") %in% names(external)))
    stop("Неожиданный формат входных таблиц; проверьте версию проекта.")

  pathology <- scores$periodontal_status == "Documented pathology"
  normal_ref <- external$sample_group == "normal_reference"
  ging_ref <- external$sample_group == "gingivitis_reference"
  if (anyNA(pathology) || length(pathology) != 25L ||
      sum(pathology) != 3L || sum(normal_ref) != 3L || sum(ging_ref) != 1L ||
      any(!is.finite(profiles$wavenumber)) ||
      any(!is.finite(profiles$contrast_all)) ||
      any(!is.finite(profiles$contrast_age)) ||
      any(!is.finite(scores$EMSC_DPT4_CH_PC1)) ||
      any(!is.finite(external$CH_PC1_pediatric_SD)))
    stop("Состав выборки или численные данные не прошли контроль.")
  if (any(diff(profiles$wavenumber) >= 0))
    stop("Ожидались убывающие волновые числа.")

  # Разметка коллег: никакой обработки сигнала при построении рисунка.
  ch2_range <- c(2915, 2935)
  marker_red <- "#CD3024"
  # Тот же увеличенный кегль заголовков, что и на рисунке 1.
  panel_title_cex <- 1.12
  if (min(profiles$wavenumber) > min(ch2_range) ||
      max(profiles$wavenumber) < max(ch2_range))
    stop("Выделяемая область C–H лежит вне спектральной оси.")

  plot_figure <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    layout(matrix(c(1, 2), nrow = 1), widths = c(1.05, 1))
    par(family = "sans", cex = 0.92)

    # А. Клинический контраст после EMSC DPT4.
    par(mar = c(4.5, 5.0, 4.8, 1.0))
    yrange <- range(c(profiles$contrast_all, profiles$contrast_age))
    ylim <- yrange + diff(yrange) * c(-0.10, 0.20)
    plot(profiles$wavenumber, profiles$contrast_all,
         type = "n", xlim = rev(range(profiles$wavenumber)), ylim = ylim,
         xlab = "Волновое число, см⁻¹",
         ylab = "Разность медианных спектров, усл. ед.", main = "",
         cex.axis = 0.94, cex.lab = 0.94)
    title(main = "А. Разность медианных спектров (C–H)", line = 3.0,
          cex.main = panel_title_cex / par("cex"))

    # Отмечаем 2915–2935 см^-1 ровно так, как показано коллегами.
    u <- par("usr")
    rect(ch2_range[1], u[3], ch2_range[2], u[4],
         col = adjustcolor("#E8C98D", alpha.f = 0.27), border = NA)
    # Все вспомогательные пунктирные линии удалены.
    # Спектральный участок 2915–2935 см⁻¹ обозначен цветовой полосой
    # и скобкой, поэтому центральная линия ему не нужна.
    lines(profiles$wavenumber, profiles$contrast_all,
          col = "#245C84", lwd = 3)
    lines(profiles$wavenumber, profiles$contrast_age,
          col = "#BD513B", lwd = 2.5, lty = 2)

    # Скобка из НАСТОЯЩИХ графических примитивов R (не ручной монтаж).
    yspan <- u[4] - u[3]
    ybr <- u[4] + 0.014 * yspan
    for (xx in ch2_range) {
      segments(xx, ybr - 0.005 * yspan, xx, ybr + 0.005 * yspan,
               col = marker_red, lwd = 1.5, xpd = NA)
    }
    segments(ch2_range[1], ybr, ch2_range[2], ybr,
             col = marker_red, lwd = 1.5, xpd = NA)
    text(mean(ch2_range), ybr + 0.033 * yspan,
         labels = expression(nu[as](CH[2])),
         col = marker_red, cex = 0.83, xpd = NA)
    legend("topleft",
           legend = c("3 против 22", "3 против 6 (возрастной контроль)"),
           col = c("#245C84", "#BD513B"), lwd = c(3, 2.5),
           lty = c(1, 2), cex = 0.8, bty = "n")
    box()

    # Б. Внешние N/N-1/N-2 и G в фиксированной PCA детей (без refit).
    par(mar = c(4.5, 9.2, 4.8, 1.0))
    child <- scores$EMSC_DPT4_CH_PC1
    refs <- external$CH_PC1_pediatric_SD
    xlim <- range(c(child, refs)) + c(-0.42, 0.52)
    plot(NA_real_, NA_real_, xlim = xlim, ylim = c(0.55, 4.45),
         xlab = "Координата C–H-PC1, SD детей", ylab = "", yaxt = "n",
         main = "", cex.axis = 0.94, cex.lab = 0.94)
    title(main = "Б. Координаты образцов по PC1",
          line = 3.0, cex.main = panel_title_cex / par("cex"))
    axis(2, at = 4:1, las = 1, cex.axis = 0.79,
         labels = c("Без записи о патологии (22)",
                    "Патология пародонта (3)",
                    "Референсы: норма (3)",
                    "Референс: гингивит (1)"))
    # Линии нуля и пунктирные разделители категорий не рисуем:
    # категории различаются подписями и формой/цветом точек.
    set.seed(20261118L) # фиксированное зерно для визуального jitter
    points(child[!pathology],
           4 + runif(sum(!pathology), -0.12, 0.12), pch = 16L,
           col = adjustcolor("#6C7A89", alpha.f = 0.78), cex = 1.05)
    points(child[pathology],
           3 + runif(sum(pathology), -0.12, 0.12), pch = 16L,
           col = "#BF4C3A", cex = 1.3)
    points(refs[normal_ref],
           2 + runif(sum(normal_ref), -0.10, 0.10), pch = 18L,
           col = "#2A729A", cex = 1.4)
    points(refs[ging_ref], 1, pch = 17L, col = "#A34E8D", cex = 1.5)
    box()
    invisible(NULL)
  }

  png_path <- file.path(results_dir,
    "Figure_CH_contrast_and_external_references_RU_Paul.png")
  pdf_path <- file.path(results_dir,
    "Figure_CH_contrast_and_external_references_RU_Paul.pdf")

  png(png_path, width = 2200, height = 1060, res = 200)
  tryCatch(plot_figure(), finally = dev.off())
  # На macOS Cairo иногда объявлен доступным, хотя библиотеки XQuartz
  # отсутствуют. В этом случае PNG готов, а PDF безопасно пропускаем.
  cairo_ready <- isTRUE(capabilities("cairo")) && (
    Sys.info()[["sysname"]] != "Darwin" ||
      all(file.exists(c("/opt/X11/lib/libSM.6.dylib",
                        "/opt/X11/lib/libXrender.1.dylib")))
  )
  if (cairo_ready) {
    grDevices::cairo_pdf(pdf_path, width = 11, height = 5.3,
                         family = "sans", onefile = TRUE)
    tryCatch(plot_figure(), finally = dev.off())
  } else {
    pdf_path <- NA_character_
    cat("Cairo недоступен: PNG создан, PDF пропущен.\n")
  }

  cat("Рисунок 2 с разметкой Павла: ", png_path, "\n", sep = "")
  if (!is.na(pdf_path)) cat("PDF: ", pdf_path, "\n", sep = "")
  invisible(c(png = png_path, pdf = pdf_path))
}

figure_output <- make_bio_figure_with_refs()
