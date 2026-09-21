# Shared preamble for the analysis stage scripts in scripts/analysis/.
# Each stage is run as `Rscript scripts/analysis/<stage>.R <analysis-name>`,
# or all of them with `Rscript scripts/run-analysis.R <analysis-name>`.
suppressPackageStartupMessages({
  library(tidyverse)
  library(logger)
  library(distionary)
  devtools::load_all(here::here(), quiet = TRUE)
})

ANALYSIS <- analysis_from_args()
cfg <- read_analysis(ANALYSIS)
out_path <- function(...) analysis_path(ANALYSIS, "out", ...)
fig_path <- function(...) analysis_path(ANALYSIS, "figs", ...)
dir.create(out_path(), showWarnings = FALSE, recursive = TRUE)
dir.create(fig_path(), showWarnings = FALSE, recursive = TRUE)

N_CORES <- max(1L, min(4L, parallel::detectCores() - 1L))

# Apply f to each cell's rows of `data` in parallel; returns a list named by
# cell_id.
map_cells <- function(data, f) {
  groups <- split(data, data$cell_id)
  res <- parallel::mclapply(groups, f, mc.cores = N_CORES, mc.preschedule = FALSE)
  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop("Cell(s) failed: ", paste(names(res)[failed], collapse = ", "), "\n",
      as.character(res[failed][[1]]), call. = FALSE)
  }
  res
}

log_info("[{ANALYSIS}] {cfg$title %||% ''}")
