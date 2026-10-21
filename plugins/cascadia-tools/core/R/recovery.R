# Simulation-recovery and power analysis for any ecomodule.
#
# run_recovery(): simulate nsim datasets from known parameters, fit each one,
# and compare estimates with the truth: bias, relative bias, RMSE, 95%
# interval coverage, interval width, convergence rate and power (the share
# of replicates whose interval excludes the null value).
#
# power_analysis(): run_recovery over a grid of designs or effect sizes and
# return power curves.
#
# Replicates run in parallel with forked workers (parallel::mclapply). Each
# replicate gets seed = base_seed + rep, so any single replicate can be
# rerun exactly with module$simulate(params, design, seed).

run_recovery <- function(module, params = NULL, design = NULL, nsim = 50,
                         settings = NULL, cores = 1, seed = 1,
                         outdir = NULL, label = NULL, verbose = TRUE) {
  params <- merge_list(module$defaults$params, params)
  design <- merge_list(module$defaults$design, design)
  settings <- merge_list(module$fast_settings, settings)
  one <- function(r) {
    s <- seed * 100000 + r
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch({
      sim <- module$simulate(params, design, s)
      f <- module$fit(sim$data, settings)
      est <- check_estimates(f$estimates)
      tr <- sim$truth
      est <- est[est$param %in% names(tr), c("param", "estimate", "sd", "lower", "upper"), drop = FALSE]
      if (!nrow(est)) stop("no estimates match names(truth): ", paste(names(tr), collapse = ","))
      est$truth <- unname(tr[est$param])
      est$rep <- r; est$seed <- s
      est$converged <- isTRUE(f$diagnostics$converged) || is.null(f$diagnostics$converged)
      est$rhat_max <- if (is.null(f$diagnostics$rhat_max)) NA_real_ else f$diagnostics$rhat_max
      est$time <- proc.time()[["elapsed"]] - t0
      est$error <- NA_character_
      est
    }, error = function(e) data.frame(param = NA_character_, estimate = NA_real_, sd = NA_real_,
                                      lower = NA_real_, upper = NA_real_, truth = NA_real_,
                                      rep = r, seed = s, converged = FALSE, rhat_max = NA_real_,
                                      time = proc.time()[["elapsed"]] - t0,
                                      error = conditionMessage(e)))
    if (verbose) message(sprintf("[%s] rep %d done (%.1fs)%s", module$id, r, out$time[1],
                                 if (!is.na(out$error[1])) paste0(" ERROR: ", out$error[1]) else ""))
    out
  }
  reps <- if (cores > 1) parallel::mclapply(seq_len(nsim), one, mc.cores = cores,
                                            mc.preschedule = TRUE)  # one fork per core, so compiled models are reused across reps
          else lapply(seq_len(nsim), one)
  bad <- vapply(reps, function(x) !is.data.frame(x), logical(1))
  if (any(bad)) reps[bad] <- lapply(which(bad), function(r)
    data.frame(param = NA_character_, estimate = NA_real_, sd = NA_real_, lower = NA_real_,
               upper = NA_real_, truth = NA_real_, rep = r, seed = seed * 100000 + r,
               converged = FALSE, rhat_max = NA_real_, time = NA_real_,
               error = "worker crashed"))
  reps <- do.call(rbind, reps)
  summ <- summarise_recovery(reps, module)
  res <- list(module = module$id, params = params, design = design, settings = settings,
              nsim = nsim, seed = seed, reps = reps, summary = summ)
  class(res) <- "eco_recovery"
  if (!is.null(outdir)) save_recovery(res, outdir, label)
  res
}

summarise_recovery <- function(reps, module) {
  ok <- reps[!is.na(reps$param), ]
  nulls <- module$null_values
  if (is.null(nulls)) nulls <- stats::setNames(rep(0, length(module$power_params)), module$power_params)
  n_fail <- length(unique(reps$rep[!is.na(reps$error)]))
  if (!nrow(ok)) return(data.frame(param = character(), n_failed = n_fail))
  sp <- split(ok, ok$param)
  out <- do.call(rbind, lapply(sp, function(d) {
    conv <- d$converged %in% TRUE
    d2 <- if (any(conv)) d[conv, ] else d
    err <- d2$estimate - d2$truth
    nv <- if (d2$param[1] %in% names(nulls)) nulls[[d2$param[1]]] else NA_real_
    data.frame(
      # n_ok: replicates summarised (the converged ones when any converged);
      # n_ran: replicates that returned estimates, converged or not.
      param = d2$param[1], truth = mean(d2$truth), n_ok = nrow(d2), n_ran = nrow(d),
      n_converged = sum(conv), mean_est = mean(d2$estimate),
      bias = mean(err), rel_bias = if (all(d2$truth != 0)) mean(err / abs(d2$truth)) else NA_real_,
      bias_mcse = stats::sd(err) / sqrt(nrow(d2)),
      rmse = sqrt(mean(err^2)),
      coverage = mean(d2$lower <= d2$truth & d2$upper >= d2$truth),
      ci_width = mean(d2$upper - d2$lower),
      mean_sd = mean(d2$sd, na.rm = TRUE), emp_sd = stats::sd(d2$estimate),
      null_value = nv,
      power = if (is.na(nv)) NA_real_ else mean(d2$lower > nv | d2$upper < nv),
      row.names = NULL)
  }))
  out$n_failed <- n_fail
  out
}

# Pass/fail gate on recovery (see protocols/acceptance-gate.md). Thresholds
# are fixed here, before any module is fit, and are not tuned per module.
# With nsim replicates the coverage MC standard error is ~sqrt(.95*.05/nsim),
# so the coverage bound is 0.95 minus 3 MC SEs, floored at 0.80.
recovery_gate <- function(summ, max_abs_rel_bias = 0.10, min_coverage = NULL,
                          min_converged = 0.9) {
  n <- summ$n_ok
  cov_lo <- if (is.null(min_coverage)) pmax(0.80, 0.95 - 3 * sqrt(0.95 * 0.05 / pmax(n, 1))) else min_coverage
  # Bias is judged against the larger of the relative tolerance and 3 MC SEs
  tol <- pmax(max_abs_rel_bias * abs(summ$truth), 3 * summ$bias_mcse)
  tol[summ$truth == 0] <- pmax(0.1 * summ$emp_sd, 3 * summ$bias_mcse)[summ$truth == 0]
  bias_ok <- abs(summ$bias) <= tol
  cov_ok <- summ$coverage >= cov_lo
  # Denominator = every replicate attempted (ran, converged or not, plus
  # crashed). Older summaries without n_ran fall back to n_ok.
  n_ran <- if (!is.null(summ$n_ran)) summ$n_ran else summ$n_ok
  conv_ok <- summ$n_converged / pmax(n_ran + summ$n_failed, 1) >= min_converged
  data.frame(param = summ$param, bias_ok = bias_ok, coverage_ok = cov_ok,
             converged_ok = conv_ok, coverage_bound = round(cov_lo, 3),
             pass = bias_ok & cov_ok & conv_ok)
}

power_analysis <- function(module, vary, what = c("design", "params"), nsim = 100,
                           params = NULL, design = NULL, settings = NULL,
                           cores = 1, seed = 1, outdir = NULL, label = "power") {
  what <- match.arg(what)
  stopifnot(length(vary) == 1)
  key <- names(vary); vals <- vary[[1]]
  rows <- lapply(seq_along(vals), function(i) {
    p <- params; d <- design
    if (what == "design") d <- merge_list(d, stats::setNames(list(vals[i]), key))
    else p <- merge_list(p, stats::setNames(list(vals[i]), key))
    r <- run_recovery(module, p, d, nsim = nsim, settings = settings, cores = cores,
                      seed = seed + i, verbose = FALSE)
    s <- r$summary; s$vary <- key; s$value <- vals[i]; s
  })
  out <- do.call(rbind, rows)
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, file.path(outdir, paste0(module$id, "_", label, ".csv")), row.names = FALSE)
  }
  out
}

save_recovery <- function(res, outdir, label = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(outdir, paste0(res$module, if (!is.null(label)) paste0("_", label)))
  utils::write.csv(res$reps, paste0(stem, "_reps.csv"), row.names = FALSE)
  utils::write.csv(res$summary, paste0(stem, "_summary.csv"), row.names = FALSE)
  saveRDS(res, paste0(stem, "_recovery.rds"))
  invisible(stem)
}

print.eco_recovery <- function(x, ...) {
  cat(sprintf("Recovery: %s, %d replicates (%d failed)\n", x$module, x$nsim,
              if (nrow(x$summary)) x$summary$n_failed[1] else x$nsim))
  s <- x$summary
  if (nrow(s)) {
    g <- recovery_gate(s)
    s$gate <- ifelse(g$pass, "PASS", "FAIL")
    num <- vapply(s, is.numeric, logical(1))
    s[num] <- lapply(s[num], signif, 3)
    print(s[, intersect(c("param", "truth", "mean_est", "bias", "rel_bias", "rmse", "coverage",
                          "ci_width", "power", "n_converged", "gate"), names(s))], row.names = FALSE)
  }
  invisible(x)
}
