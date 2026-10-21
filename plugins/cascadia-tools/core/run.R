# CascadiaTools command line.
#
#   scripts/Rscript core/run.R list
#   scripts/Rscript core/run.R describe <module>
#   scripts/Rscript core/run.R selftest <module> [--nsim 10] [--cores 4]
#   scripts/Rscript core/run.R simulate <module> [--config cfg.yaml] [--seed 1] [--out dir]
#   scripts/Rscript core/run.R fit      <module> --config cfg.yaml [--out dir]
#   scripts/Rscript core/run.R figures  <module> [--fit dir/<module>_fit.rds] [--config cfg.yaml] [--out dir]
#   scripts/Rscript core/run.R recover  <module> [--config cfg.yaml] [--nsim 100] [--cores 8] [--out dir]
#   scripts/Rscript core/run.R power    <module> --vary n_sites=50,100,200 [--what design|params]
#                               [--config cfg.yaml] [--nsim 100] [--cores 8] [--out dir]
#
# A config YAML may hold: params, design, settings (overrides of the module
# defaults) and data (module-specific paths/column names, used by `fit`).
# Outputs go to --out (default: results/<module>/ under the current directory).

local({
  root <- Sys.getenv("CASCADIATOOLS_ROOT", Sys.getenv("ECOTOOLS_ROOT", ""))
  if (!nzchar(root)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
    root <- if (!is.na(f)) dirname(dirname(normalizePath(f))) else getwd()
  }
  Sys.setenv(CASCADIATOOLS_ROOT = normalizePath(root, mustWork = FALSE))
})
source(file.path(Sys.getenv("CASCADIATOOLS_ROOT"), "core", "load.R"))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) { cat(readLines(file.path(eco_root(), "core", "run.R"))[3:13], sep = "\n"); quit(status = 1) }
cmd <- args[1]
pos <- args[-1][!grepl("^--", args[-1]) & !(seq_along(args[-1]) %in% (which(grepl("^--", args[-1])) + 1))]
opt <- function(name, default = NULL) {
  i <- which(args == paste0("--", name))
  if (length(i)) args[i[1] + 1] else default
}
num <- function(name, default) as.numeric(opt(name, default))

read_config <- function() {
  f <- opt("config")
  if (is.null(f)) return(list())
  eco_require("yaml")
  cfg <- yaml::read_yaml(f)
  cfg$.path <- normalizePath(f)
  cfg$.dir <- dirname(cfg$.path)
  cfg
}

if (cmd == "list") {
  for (id in list_modules()) {
    m <- tryCatch(load_module(id), error = function(e) NULL)
    cat(sprintf("%-22s %s\n", id, if (is.null(m)) "(failed to load)" else m$title))
  }
  quit(status = 0)
}

id <- pos[1]
if (is.na(id)) stop("module id required")
m <- load_module(id)
cfg <- read_config()
out <- opt("out", file.path(getwd(), "results", id))  # results go where the researcher is working, never into the toolkit
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cores <- num("cores", max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))))

if (cmd == "describe") {
  cat(m$title, "\n\n", m$description, "\n\nEngine: ", m$engine, "\n\nDefault params:\n", sep = "")
  str(m$defaults$params); cat("\nDefault design:\n"); str(m$defaults$design)
  cat("\nPower parameters:", m$power_params, "\n")
  if (length(m$references)) cat("\nReferences:\n", paste("-", m$references), sep = "\n")
} else if (cmd == "simulate") {
  p <- merge_list(m$defaults$params, cfg$params); d <- merge_list(m$defaults$design, cfg$design)
  sim <- m$simulate(p, d, num("seed", 1))
  f <- file.path(out, sprintf("%s_sim_seed%d.rds", id, num("seed", 1)))
  saveRDS(sim, f)
  if (is.function(m$write_data)) m$write_data(sim$data, out)
  cat("Simulated data written to", f, "\nTruth:\n"); print(sim$truth)
} else if (cmd == "fit") {
  if (!is.function(m$read_data)) stop("module ", id, " has no read_data(); fit simulated data with `selftest`")
  data <- m$read_data(cfg)
  if (is.function(m$check_data)) {
    pr <- m$check_data(data)
    if (length(pr)) { cat("DATA PROBLEMS:\n", paste("-", pr), sep = "\n"); quit(status = 2) }
  }
  st <- merge_list(m$settings, cfg$settings)
  f <- m$fit(data, st)
  print(f$estimates, row.names = FALSE)
  utils::write.csv(f$estimates, file.path(out, paste0(id, "_estimates.csv")), row.names = FALSE)
  saveRDS(f, file.path(out, paste0(id, "_fit.rds")))
  saveRDS(list(data = data, settings = st), file.path(out, paste0(id, "_data.rds")))
  rep <- analysis_report(m, f, data, st, file.path(out, paste0(id, "_report.html")), opt("config"))
  cat("Report:", rep, "\n")
  fx <- tryCatch(export_figures(m, f, data, fig_opts(cfg$figures), file.path(out, "figures")),
                 error = function(e) { message("figure export failed: ", conditionMessage(e)); NULL })
  if (!is.null(fx)) cat("Figures:", file.path(out, "figures", "gallery.html"), "\n")
} else if (cmd == "figures") {
  # Re-draw publication figures from a saved fit without refitting. Edit the
  # `figures:` block of the config (see skills/eco-figures) and re-run.
  ff <- opt("fit", file.path(out, paste0(id, "_fit.rds")))
  f <- readRDS(ff)
  dd <- readRDS(sub("_fit\\.rds$", "_data.rds", ff))
  fx <- export_figures(m, f, dd$data, fig_opts(cfg$figures), opt("figdir", file.path(dirname(ff), "figures")))
  print(fx, row.names = FALSE)
  cat("Gallery:", file.path(opt("figdir", file.path(dirname(ff), "figures")), "gallery.html"), "\n")
} else if (cmd %in% c("recover", "selftest")) {
  nsim <- num("nsim", if (cmd == "selftest") 4 else 100)
  r <- run_recovery(m, cfg$params, cfg$design, nsim = nsim, settings = cfg$settings,
                    cores = cores, seed = num("seed", 1), outdir = out, label = cmd)
  print(r)
  if (nrow(r$summary)) print(recovery_gate(r$summary), row.names = FALSE)
  rep <- recovery_report(m, r, file.path(out, paste0(id, "_", cmd, "_report.html")))
  cat("Report:", rep, "\n")
} else if (cmd == "power") {
  v <- opt("vary"); if (is.null(v)) stop("--vary key=v1,v2,... required")
  kv <- strsplit(v, "=")[[1]]
  vals <- as.numeric(strsplit(kv[2], ",")[[1]])
  pw <- power_analysis(m, stats::setNames(list(vals), kv[1]), what = opt("what", "design"),
                       nsim = num("nsim", 100), params = cfg$params, design = cfg$design,
                       settings = cfg$settings, cores = cores, seed = num("seed", 1), outdir = out)
  print(pw[, c("vary", "value", "param", "power", "coverage", "rel_bias")], row.names = FALSE)
  rep <- power_report(m, pw, file.path(out, paste0(id, "_power_report.html")),
                      params = merge_list(m$defaults$params, cfg$params),
                      design = merge_list(m$defaults$design, cfg$design), nsim = num("nsim", 100))
  cat("Report:", rep, "\n")
} else stop("unknown command: ", cmd)
