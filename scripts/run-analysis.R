# Run an analysis end to end, or selected stages.
#
#   Rscript scripts/run-analysis.R <analysis-name>            # all stages
#   Rscript scripts/run-analysis.R <analysis-name> 2 4 5      # stages 2, 4, 5
#
# Stages (scripts/analysis/): 1 training, 2 fit, 3 marginals, 4 queries,
# 5 figures. Requires the shared data from scripts 1-3 in derived/.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/run-analysis.R <analysis-name> [stages]")
name <- args[[1]]
stage_files <- sort(list.files(here::here("scripts", "analysis"), pattern = "^[0-9]-.*\\.R$"))
wanted <- if (length(args) > 1) as.integer(args[-1]) else seq_along(stage_files)
rscript <- file.path(R.home("bin"), "Rscript")
for (f in stage_files[as.integer(substr(stage_files, 1, 1)) %in% wanted]) {
  message("\n=== ", name, ": ", f, " ===")
  t0 <- Sys.time()
  status <- system2(rscript, c(here::here("scripts", "analysis", f), name))
  if (status != 0) stop("Stage ", f, " failed for ", name, call. = FALSE)
  message(sprintf("=== %s finished in %.1f min", f, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
