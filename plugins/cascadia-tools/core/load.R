# Source this to use CascadiaTools interactively:
#   source("<path to CascadiaTools>/core/load.R")
#   m <- load_module("occu-single"); r <- run_recovery(m, nsim = 20, cores = 4)
local({
  here <- Sys.getenv("CASCADIATOOLS_ROOT", Sys.getenv("ECOTOOLS_ROOT", ""))
  if (!nzchar(here)) {
    # When sourced, source() records this file's path as `ofile` in a frame.
    of <- NULL
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) { of <- f; break }
    }
    here <- if (!is.null(of)) dirname(dirname(normalizePath(of))) else getwd()
  }
  Sys.setenv(CASCADIATOOLS_ROOT = normalizePath(here, mustWork = FALSE))
  for (f in c("env.R", "contract.R", "mcmc.R", "recovery.R", "report.R", "figures.R", "glm_helpers.R"))
    if (file.exists(file.path(here, "core", "R", f))) sys.source(file.path(here, "core", "R", f), envir = globalenv())
})
