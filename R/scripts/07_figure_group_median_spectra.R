# Рисунок 1. Редакция Павла: спектроскопические скобки заданы кодом,
# а не нанесены поверх готового PNG. Данные и масштабы не изменены.
# Рисунок на РУССКОМ ЯЗЫКЕ: реальные ГРУППОВЫЕ МЕДИАНЫ EMSC DPT4
# спектров детей в двух областях: C-H (положительный поисковый результат)
# и Protein/Amide (область сравнения).
#
# ВАЖНО: никаких повторных EMSC, PCA, подгонки моделей или нового тестирования.
# Читаются готовые CSV из проверенного запуска 00_run_biology.R.
# Три случая с документированной патологией; у 22 других детей
# НЕ записана патология пародонта (они не объявляются здоровой нормой).
# Возрастной контроль из 6 детей является подмножеством этих 22.
#
# На Mac: в RStudio выполните source(file.choose()) и выберите ЭТОТ скрипт.
# Получите PNG и PDF в biology_v2/results/.

find_biology_dir_for_spectral_figure <- function() {
  explicit_root <- getOption("gcf.figure.biology_dir", NULL)
  if (!is.null(explicit_root))
    return(normalizePath(explicit_root, mustWork = TRUE))
  frames <- sys.frames()
  files <- unlist(lapply(frames, function(f) f$ofile), use.names = FALSE)
  files <- files[!is.na(files) &
                   basename(files) == "07_figure_group_median_spectra.R"]
  if (length(files)) {
    script <- normalizePath(tail(files, 1L), mustWork = TRUE)
    return(normalizePath(file.path(dirname(script), ".."), mustWork = TRUE))
  }
  # Never read an old global biology_dir left by a previous RStudio run.
  stop("Запустите скрипт через source(file.choose()) из RStudio.")
}

read_checked_group_profile <- function(path, bounds) {
  if (!file.exists(path)) stop("Нет проверенного файла результатов: ", path,
                               "\nСначала выполните 00_run_biology.R.")
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("wavenumber", "median_case", "median_control_all",
              "median_control_age", "contrast_all", "contrast_age")
  if (!all(needed %in% names(d)))
    stop("Неподходящий формат таблицы: ", basename(path))
  if (nrow(d) < 25L || any(!is.finite(as.matrix(d[, needed]))))
    stop("Неполные или нечисловые спектры: ", basename(path))
  if (any(diff(d$wavenumber) >= 0) || anyDuplicated(d$wavenumber))
    stop("Волновые числа должны убывать без повторов: ", basename(path))
  if (min(d$wavenumber) < bounds[1L] - 1 ||
      max(d$wavenumber) > bounds[2L] + 1)
    stop("В таблице не та спектральная область: ", basename(path))
  if (max(abs(d$median_case - d$median_control_all - d$contrast_all)) > 1e-8 ||
      max(abs(d$median_case - d$median_control_age - d$contrast_age)) > 1e-8)
    stop("Медианы не согласуются с проверенными контрастами: ", basename(path))
  d
}

make_group_median_spectral_figure <- function() {
  root <- find_biology_dir_for_spectral_figure()
  result_dir <- file.path(root, "results")
  if (!dir.exists(result_dir))
    stop("Нет папки results. Сначала выполните 00_run_biology.R.")

  # Только основной, количественно проверенный вариант EMSC DPT4.
  ch <- read_checked_group_profile(
    file.path(result_dir, "03_EMSC_DPT4_CH_group_profiles.csv"),
    c(2870, 3050))
  protein <- read_checked_group_profile(
    file.path(result_dir, "03_EMSC_DPT4_Protein_group_profiles.csv"),
    c(1540, 1720))

  # Проверяем, что папка действительно принадлежит нашей версии анализа,
  # а не более ранним исследованиям с другой коррекцией носителя.
  summary_path <- file.path(result_dir, "02_pca_region_summary.csv")
  if (!file.exists(summary_path)) stop("Отсутствует 02_pca_region_summary.csv.")
  s <- read.csv(summary_path, stringsAsFactors = FALSE)
  for (r in c("CH", "Protein")) {
    row <- s[s$carrier_reference == "EMSC_DPT4" & s$region == r, , drop = FALSE]
    expected <- if (r == "CH") 76.534297475 else 79.321984962
    if (nrow(row) != 1L ||
        !is.finite(row$PC1_variance_pct) ||
        abs(row$PC1_variance_pct - expected) > 0.02)
      stop("Файлы results не соответствуют проверенной EMSC DPT4: ", r)
  }

  colors <- c(case = "#BB513D", all = "#245C84", age = "#59616D")
  # Единый увеличенный кегль заголовков панелей для обоих рисунков.
  panel_title_cex <- 1.12

  # Спектроскопическая разметка ПОЛНОСТЬЮ ВОСПРОИЗВОДИМА:
  # ширины интервалов и подписи зафиксированы в коде.
  # Они обозначают области колебаний, а не точную идентификацию веществ.
  annotations <- list(
    CH2 = c(2915, 2935),
    Amide_I = c(1600, 1700),
    Amide_II = c(1540, 1580)
  )
  annotation_col <- "#CD3024"

  add_spectral_bracket <- function(limits, label) {
    u <- par("usr")
    yr <- u[4] - u[3]
    # Маленькая скобка над верхней границей области построения.
    ybar <- u[4] + 0.014 * yr
    tick <- 0.010 * yr
    # Если граница окна проходит между отсчётами FTIR, ограничиваем
    # положение скобки фактическим диапазоном панели.
    xx <- sort(pmax(min(u[1:2]), pmin(max(u[1:2]), limits)))
    segments(xx[1], ybar, xx[2], ybar,
             col = annotation_col, lwd = 1.5, xpd = NA)
    segments(xx[1], ybar - tick/2, xx[1], ybar + tick/2,
             col = annotation_col, lwd = 1.5, xpd = NA)
    segments(xx[2], ybar - tick/2, xx[2], ybar + tick/2,
             col = annotation_col, lwd = 1.5, xpd = NA)
    text(mean(xx), ybar + 0.034 * yr, labels = label,
         cex = 0.78, col = annotation_col, xpd = NA)
  }

  draw_region <- function(d, heading, highlight = FALSE, show_legend = FALSE) {
    ys <- c(d$median_case, d$median_control_all, d$median_control_age)
    span <- diff(range(ys))
    padding <- max(0.08 * span, 1e-8)
    ylim <- range(ys) + c(-padding, padding)
    plot(d$wavenumber, d$median_case, type = "n",
         xlim = rev(range(d$wavenumber)), ylim = ylim,
         main = "",
         xlab = "Волновое число, см⁻¹",
         ylab = "Сигнал после EMSC, усл. ед.",
         cex.lab = 0.91, cex.axis = 0.84,
         xaxs = "i")
    title(main = heading, line = 3.05,
          cex.main = panel_title_cex / par("cex"))
    # Вспомогательные пунктирные линии сетки и нулевой уровень
    # не отображаем; пунктирную КРИВУЮ возрастного контроля сохраняем.
    if (highlight) {
      # Отмечаем регион интерпретации, НЕ меняя масштаба y или интенсивностей.
      rect(2915, par("usr")[3], 2935, par("usr")[4],
           col = adjustcolor("#E8C98D", alpha.f = 0.27), border = NA)
    }
    # Нулевая горизонтальная центрирующая линия намеренно не отображается.
    lines(d$wavenumber, d$median_control_all,
          col = colors["all"], lwd = 2.7)
    lines(d$wavenumber, d$median_control_age,
          col = colors["age"], lwd = 2.4, lty = 2)
    lines(d$wavenumber, d$median_case,
          col = colors["case"], lwd = 3.0)
    if (highlight) {
      add_spectral_bracket(annotations$CH2,
                            expression(nu[as](CH[2])))
    } else {
      add_spectral_bracket(annotations$Amide_I, "Амид I")
      add_spectral_bracket(annotations$Amide_II, "Амид II")
    }
    if (show_legend) {
      legend("topleft",
             legend = c("Патология пародонта (n = 3)",
                        "Без записи о патологии (n = 22)",
                        "Возрастной контроль (n = 6)"),
             col = colors[c("case", "all", "age")],
             lty = c(1, 1, 2), lwd = c(3, 2.7, 2.4),
             cex = 0.74, bty = "n", seg.len = 1.8)
    }
    box()
  }

  draw_figure <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(1, 2), mar = c(4.7, 5.2, 5.0, 1.05),
        oma = c(0, 0, 0, 0), family = "sans", cex = 0.95)
    draw_region(ch, "А. Область C–H: 2870–3050 см⁻¹", highlight = TRUE,
                show_legend = TRUE)
    draw_region(protein, "Б. Амидная область: 1540–1720 см⁻¹")
  }

  # Не перезаписываем прежние, ещё не размеченные варианты.
  png_path <- file.path(result_dir, "Figure_group_median_spectra_CH_and_Amide_RU_Paul.png")
  pdf_path <- file.path(result_dir, "Figure_group_median_spectra_CH_and_Amide_RU_Paul.pdf")

  png(png_path, width = 2300, height = 1100, res = 230)
  tryCatch(draw_figure(), finally = dev.off())
  # Встроенный pdf() с Helvetica не гарантирует кириллицу; используем Cairo.
  # На macOS без библиотек XQuartz пропускаем только необязательный PDF.
  cairo_ready <- isTRUE(capabilities("cairo")) && (
    Sys.info()[["sysname"]] != "Darwin" ||
      all(file.exists(c("/opt/X11/lib/libSM.6.dylib",
                        "/opt/X11/lib/libXrender.1.dylib")))
  )
  if (cairo_ready) {
    grDevices::cairo_pdf(pdf_path, width = 10, height = 4.8,
                         family = "sans", onefile = TRUE)
    tryCatch(draw_figure(), finally = dev.off())
  } else {
    cat("Cairo недоступен: русский PDF не создан, PNG полностью готов.\n")
    pdf_path <- NA_character_
  }

  cat("Готово: ДВЕ панели с медианными EMSC DPT4 спектрами.\n")
  cat("PNG: ", png_path, "\n", sep = "")
  if (!is.na(pdf_path)) cat("PDF: ", pdf_path, "\n", sep = "")
  cat("Панели имеют НЕЗАВИСИМЫЕ шкалы y; внутри каждой панели ",
      "все группы показаны в ОДИНАКОВОМ масштабе.\n", sep = "")
  invisible(c(png = png_path, pdf = pdf_path))
}

figure_files <- make_group_median_spectral_figure()
