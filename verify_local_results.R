# Verify numeric checkpoints and final figures from the current locally authorized run.
find_verify_root <- function() {
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names=FALSE)
  paths <- paths[!is.na(paths) & basename(paths) == "verify_local_results.R"]
  if (length(paths)) return(dirname(normalizePath(tail(paths, 1L), mustWork=TRUE)))
  args <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (length(args)) return(dirname(normalizePath(sub("^--file=", "", tail(args, 1L)), mustWork=TRUE)))
  stop("Запустите через source(file.choose()) или Rscript.")
}
root_dir <- find_verify_root()
r <- file.path(root_dir, "R", "results")
p <- file.path(r, "05_independent_numeric_checks.csv")
if (!file.exists(p)) stop("Нет 05_independent_numeric_checks.csv. Сначала запустите run_all.R.")
tab <- read.csv(p, stringsAsFactors = FALSE)
if (!all(c("metric","independent_reference","observed_R","absolute_tolerance","pass") %in% names(tab)) || nrow(tab) != 5L)
  stop("Неожиданный формат численных проверок.")
if (!all(tab$pass) || anyNA(tab$pass)) stop("Есть несовпавшие контрольные значения: ", p)
if (any(abs(tab$independent_reference-tab$observed_R) > tab$absolute_tolerance))
  stop("Данные проверок не согласуются с допусками.")
for (name in c("Figure_group_median_spectra_CH_and_Amide_RU_Paul.png", "Figure_CH_contrast_and_external_references_RU_Paul.png")) {
  q <- file.path(r, name)
  if (!file.exists(q) || file.info(q)$size < 10000)
    stop("Отсутствует или поврежден один из окончательных рисунков: ", name)
}
m <- file.path(r, "04_NMF_master_summary.csv")
if (!file.exists(m)) stop("Нет summary новой NMF: ", m)
cat("PASS: 5/5 numeric checks, NMF summary present, 2 final PNG present in ", r, "\n", sep = "")
