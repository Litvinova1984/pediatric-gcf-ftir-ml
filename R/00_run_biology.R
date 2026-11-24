# Биологически ориентированная ветка статьи «Сорбционные процессы».
# Один запуск из чистой R-сессии, без внешних R-пакетов.
# Вход: data/public/GCF_MASTER_34_PUBLIC.csv (предпочтительный) или data/private/GCF_MASTER_34.csv
# Выход: biology_v2/results/ (не заменяет старые результаты).

identify_script_dir <- function() {
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  paths <- paths[!is.na(paths) & basename(paths) == "00_run_biology.R"]
  if (length(paths)) return(dirname(normalizePath(tail(paths, 1L), mustWork = TRUE)))
  file_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_args)) return(dirname(normalizePath(sub("^--file=", "", tail(file_args, 1L)), mustWork = TRUE)))
  stop("Запускайте через source(абсолютный_путь_к_00_run_biology.R) или Rscript.")
}

run_biology <- function() {
  biology_dir <- identify_script_dir()
  project_dir <- normalizePath(file.path(biology_dir, ".."), mustWork = TRUE)
  public_file <- file.path(project_dir, "data", "public", "GCF_MASTER_34_PUBLIC.csv")
  private_file <- file.path(project_dir, "data", "private", "GCF_MASTER_34.csv")
  if (file.exists(public_file)) {
    input_file <- public_file
  } else if (file.exists(private_file)) {
    input_file <- private_file
  } else {
    stop("Входная матрица не найдена: ожидается data/public/GCF_MASTER_34_PUBLIC.csv или data/private/GCF_MASTER_34.csv.")
  }
  result_dir <- file.path(biology_dir, "results")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  logcon <- file(file.path(result_dir, "RUN_BIOLOGY_LOG.txt"), open = "wt")
  previous_sinks <- sink.number()
  sink(logcon, split = TRUE)
  on.exit({
    while (sink.number() > previous_sinks) sink()
    close(logcon)
  }, add = TRUE)
  cat("Новая биологическая ветка GCF, запуск:", format(Sys.time()), "\n")
  cat("Отдельная версия: прежние результаты и скрипты не перезаписываются.\n")
  script_names <- c("01_inputs_emsc.R", "02_pca_regions.R",
                    "03_biological_contrasts.R", "04_nmf_master.R",
                    "05_write_results.R")
  skip_nmf <- identical(Sys.getenv("GCF_BIO_SKIP_NMF"), "1")
  for (script in script_names) {
    if (skip_nmf && script == "04_nmf_master.R") {
      cat("NMF пропущена по GCF_BIO_SKIP_NMF=1.\n")
      next
    }
    cat("\nЭтап:", script, "\n")
    source(file.path(biology_dir, "scripts", script), local = TRUE,
           chdir = FALSE, encoding = "UTF-8")
  }
  cat("\nГОТОВО. Результаты:", result_dir, "\n")
  invisible(result_dir)
}

run_biology()
