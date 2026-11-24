# Generate both publication figures strictly from this repository's R/results.
# Can be launched standalone via source(file.choose()) or by root run_all.R.

find_paul_script_dir <- function() {
  explicitly_set <- getOption("gcf.figure.script_dir", NULL)
  if (!is.null(explicitly_set))
    return(normalizePath(explicitly_set, mustWork = TRUE))
  paths <- unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE)
  paths <- paths[!is.na(paths) & basename(paths) == "00_make_Paul_figures.R"]
  if (length(paths))
    return(dirname(normalizePath(tail(paths, 1L), mustWork = TRUE)))
  stop("Запустите 00_make_Paul_figures.R через source(file.choose()).")
}

make_both_paul_figures <- function() {
  script_dir <- find_paul_script_dir()
  target_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
  expected_dir <- getOption("gcf.figure.biology_dir", NULL)
  if (!is.null(expected_dir) &&
      !identical(normalizePath(expected_dir, mustWork = TRUE), target_dir))
    stop("Несогласованная папка рисунков: каталог скриптов не относится к текущему проекту.")
  results_dir <- file.path(target_dir, "results")
  if (!dir.exists(results_dir))
    stop("Нет папки R/results. Сначала запустите run_all.R или R/00_run_biology.R.")
  previous <- getOption("gcf.figure.biology_dir", NULL)
  on.exit(options(gcf.figure.biology_dir = previous), add = TRUE)
  options(gcf.figure.biology_dir = target_dir)

  source(file.path(script_dir, "06_figure_with_references.R"),
         local = TRUE, encoding = "UTF-8")
  source(file.path(script_dir, "07_figure_group_median_spectra.R"),
         local = TRUE, encoding = "UTF-8")

  expected <- file.path(results_dir, c(
    "Figure_CH_contrast_and_external_references_RU_Paul.png",
    "Figure_group_median_spectra_CH_and_Amide_RU_Paul.png"))
  actual <- c(unname(figure_output["png"]), unname(figure_files["png"]))
  if (length(actual) != 2L || anyNA(actual) ||
      any(!file.exists(expected)) ||
      !identical(normalizePath(actual, mustWork = TRUE),
                 normalizePath(expected, mustWork = TRUE)))
    stop("Ошибка путей: рисунки не созданы внутри текущего R/results.")
  cat("PASS: оба рисунка построены в ", results_dir, "\n", sep = "")
}

make_both_paul_figures()
