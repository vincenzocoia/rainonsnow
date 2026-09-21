#' Locate an analysis folder and its files
#'
#' An analysis lives in `analyses/<name>/`, with its settings in
#' `analysis.yaml`, model outputs in `out/` and figures in `figs/`. Shared
#' upstream data (hourly table, POT peaks) stays in `derived/`.
#'
#' @param name Analysis name: the folder name under `analyses/`.
#' @param ... Path components appended to the analysis folder.
#' @param root Repository root.
#' @returns A file path.
#' @examples
#' \dontrun{
#' analysis_path("rain-llqr", "out", "models.rds")
#' }
#' @export
analysis_path <- function(name, ..., root = here::here()) {
  checkmate::assert_string(name)
  file.path(root, "analyses", name, ...)
}

#' List the analyses in the repository
#'
#' @param root Repository root.
#' @returns A character vector of analysis names (folders containing an
#'   `analysis.yaml`).
#' @export
list_analyses <- function(root = here::here()) {
  dirs <- list.dirs(file.path(root, "analyses"), recursive = FALSE)
  has_cfg <- file.exists(file.path(dirs, "analysis.yaml"))
  sort(basename(dirs[has_cfg]))
}

#' Read and validate an analysis configuration
#'
#' Fills in defaults so downstream stages can rely on every field being
#' present. See `analyses/README.md` for the meaning of each field.
#'
#' @param name Analysis name.
#' @param root Repository root.
#' @returns A list: the configuration, with `name` added.
#' @export
read_analysis <- function(name, root = here::here()) {
  path <- analysis_path(name, "analysis.yaml", root = root)
  if (!file.exists(path)) {
    rlang::abort(paste0(
      "No analysis called '", name, "' (looked for ", path, "). ",
      "Available: ", paste(list_analyses(root), collapse = ", ")
    ))
  }
  cfg <- yaml::read_yaml(path)
  cfg$name <- name
  cfg$response <- cfg$response %||% "runoff_hourly"
  cfg$predictors <- unlist(cfg$predictors, use.names = FALSE)
  checkmate::assert_character(cfg$predictors, min.len = 1, any.missing = FALSE)
  unknown <- setdiff(cfg$predictors, names(feature_registry()))
  if (length(unknown) > 0) {
    rlang::abort(paste0(
      "Unknown predictor(s): ", paste(unknown, collapse = ", "),
      ". Add them to feature_registry() in R/features.R."
    ))
  }
  checkmate::assert_list(cfg$model)
  checkmate::assert_choice(cfg$model$type, names(dl_model_registry()))
  cfg$model$args <- cfg$model$args %||% list()
  cfg$tail <- utils::modifyList(
    list(type = "gp", adaptive_threshold = 0.5),
    cfg$tail %||% list()
  )
  cfg$drivers <- utils::modifyList(
    list(marginal_family = "gamma", family_set = "parametric"),
    cfg$drivers %||% list()
  )
  cfg$queries <- utils::modifyList(
    list(
      return_periods = c(2, 5, 10, 20, 50, 100),
      trigger = list(probs = c(0.1, 0.5, 0.9)),
      grid = list(n = if (length(cfg$predictors) == 1) 200L else 60L, mult = 1.1)
    ),
    cfg$queries %||% list()
  )
  cfg$queries$return_periods <- unlist(cfg$queries$return_periods)
  cfg$queries$trigger$probs <- unlist(cfg$queries$trigger$probs)
  cfg$queries$trigger$target <- cfg$queries$trigger$target %||% cfg$predictors[1]
  checkmate::assert_choice(cfg$queries$trigger$target, cfg$predictors)
  cfg
}

#' Analysis name from the command line
#'
#' Stage scripts are run as `Rscript scripts/analysis/<stage>.R <name>`. When
#' sourced interactively, set `ANALYSIS <- "<name>"` in the global environment
#' first.
#'
#' @returns The analysis name.
#' @export
analysis_from_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    return(args[[1]])
  }
  if (exists("ANALYSIS", envir = globalenv())) {
    return(get("ANALYSIS", envir = globalenv()))
  }
  rlang::abort(
    "Give the analysis name: `Rscript <stage>.R <name>`, or set ANALYSIS first."
  )
}

#' Pick the focus cell for single-cell figures
#'
#' Uses `focus_cell` from the configuration when it exists in the data;
#' otherwise the cell with the most peaks where every predictor is positive
#' (for rain + snowmelt, the "mixed regime" cell).
#'
#' @param cfg Analysis configuration from [read_analysis()].
#' @param training Training table with `cell_id` and the predictor columns.
#' @returns A single cell id.
#' @export
choose_focus_cell <- function(cfg, training) {
  cells <- sort(unique(training$cell_id))
  fc <- cfg$focus_cell
  if (!is.null(fc) && fc %in% cells) {
    return(as.integer(fc))
  }
  pos <- rep(TRUE, nrow(training))
  for (p in cfg$predictors) pos <- pos & training[[p]] > 0
  counts <- table(training$cell_id[pos])
  if (length(counts) == 0) {
    return(cells[[1]])
  }
  as.integer(names(counts)[which.max(counts)])
}
