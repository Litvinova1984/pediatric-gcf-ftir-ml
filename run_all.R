# Reproduce the public analysis and two publication figures in the SAME repository.
# Usage: Rscript run_all.R OR source(file.choose(), encoding = "UTF-8")

find_release_root <- function(name) {
  frames <- sys.frames()
  paths <- unlist(lapply(frames, function(x) x$ofile), use.names = FALSE)
  paths <- paths[!is.na(paths) & basename(paths) == name]
  if (length(paths))
    return(dirname(normalizePath(tail(paths, 1L), mustWork = TRUE)))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args))
    return(dirname(normalizePath(sub("^--file=", "", tail(args, 1L)), mustWork = TRUE)))
  stop("Запустите run_all.R через source(file.choose()) или Rscript.")
}

run_gcf_public_release <- function() {
  root_dir <- find_release_root("run_all.R")
  script_dir <- file.path(root_dir, "R", "scripts")
  biology_dir <- normalizePath(file.path(root_dir, "R"), mustWork = TRUE)
  public_input <- file.path(root_dir, "data", "public", "GCF_MASTER_34_PUBLIC.csv")
  private_input <- file.path(root_dir, "data", "private", "GCF_MASTER_34.csv")
  if (!file.exists(public_input) && !file.exists(private_input))
    stop("Нет входной матрицы в data/public/ или data/private/.")

  # The scientific pipeline writes to R/results in this release.
  source(file.path(biology_dir, "00_run_biology.R"), local = TRUE, encoding = "UTF-8")

  # Explicitly lock figure paths to this release, not stale RStudio globals.
  previous_dir <- getOption("gcf.figure.biology_dir", NULL)
  previous_scripts <- getOption("gcf.figure.script_dir", NULL)
  on.exit(options(gcf.figure.biology_dir = previous_dir,
                  gcf.figure.script_dir = previous_scripts), add = TRUE)
  options(gcf.figure.biology_dir = biology_dir,
          gcf.figure.script_dir = script_dir)
  source(file.path(script_dir, "00_make_Paul_figures.R"),
         local = TRUE, encoding = "UTF-8")
  cat("ГОТОВО: анализ и оба рисунка выполнены из текущего репозитория: ",
      root_dir, "\n", sep = "")
}

run_gcf_public_release()
