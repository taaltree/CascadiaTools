# Shared helpers for the fisheries-stock module.

`%||%` <- function(a, b) if (is.null(a)) b else a

fs_fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
fs_sig <- function(x, d = 3) format(signif(x, d), big.mark = ",", scientific = FALSE, trim = TRUE)
fs_pct <- function(p) ifelse(p < 0.005, "<1%", ifelse(p > 0.995, ">99%", sprintf("%.0f%%", 100 * p)))

# Principal branch of the Lambert W function (x >= 0), Newton iterations.
# Used for the exact Ricker S_MSY / U_MSY (Scheuerell 2016).
lambert_w0 <- function(x) {
  w <- ifelse(x < 1, x * (1 - x + 1.5 * x^2), log(x) - log(pmax(log(x), 1e-12)))
  w[!is.finite(w) | w < 0] <- log1p(x[!is.finite(w) | w < 0])
  for (i in 1:60) {
    ew <- exp(w)
    f <- w * ew - x
    w <- w - f / (ew * (w + 1) - (w + 2) * f / (2 * w + 2))   # Halley step
  }
  w
}

# Summarise per-chain matrices of (derived) draws into contract rows with
# R-hat and bulk ESS. chains: list of matrices with identical columns.
fs_summarise_chains <- function(chains, prob = 0.95) {
  s <- summarise_draws_eco(chains, prob = prob)
  s$estimates
}

# Rows from a matrix of Monte-Carlo draws (no chains): parametric bootstrap /
# MVN-propagated likelihood intervals.
fs_draw_rows <- function(draws, est = NULL, prob = 0.95) {
  draws <- as.matrix(draws)
  a <- (1 - prob) / 2
  q <- function(v) if (all(is.na(v))) c(NA, NA) else stats::quantile(v, c(a, 1 - a), na.rm = TRUE, names = FALSE)
  qs <- t(apply(draws, 2, q))
  data.frame(param = colnames(draws),
             estimate = if (is.null(est)) colMeans(draws, na.rm = TRUE) else unname(est[colnames(draws)]),
             sd = apply(draws, 2, stats::sd, na.rm = TRUE),
             lower = qs[, 1], upper = qs[, 2], row.names = NULL)
}

# Standard rows-from-table helper with the extra columns present everywhere.
fs_rows <- function(param, estimate, sd, lower, upper, scale = "", rhat = NA_real_, ess_bulk = NA_real_) {
  data.frame(param = param, estimate = estimate, sd = sd, lower = lower, upper = upper,
             rhat = rhat, ess_bulk = ess_bulk, row.names = NULL)
}

fs_col <- function(df, col, what) {
  if (is.null(col) || !col %in% names(df))
    stop(sprintf("Column '%s' (%s) not found in the data file. Columns present: %s. Set data: %s_col in the config.",
                 col %||% "<unset>", what, paste(names(df), collapse = ", "), what), call. = FALSE)
  df[[col]]
}

# Convert a data column to numbers, failing with a message a non-coder can act
# on (which column, which rows, which values) instead of silently producing
# NA or a text column. Empty cells and NA stay NA.
fs_num <- function(x, col, what) {
  if (is.numeric(x)) return(as.numeric(x))
  s <- trimws(as.character(x))
  s[s %in% c("", "NA", "na", "N/A", ".")] <- NA
  v <- suppressWarnings(as.numeric(s))
  bad <- which(!is.na(s) & is.na(v))
  if (length(bad)) {
    hint <- if (any(grepl("^[0-9]+\\+$", s[bad]))) " Plus-group ages such as '10+' cannot be used in a catch curve: give the age as a number or remove those rows."
            else if (any(grepl("[0-9],[0-9]", s[bad]))) " Remove thousands separators (write 1234.5, not 1,234.5)."
            else ""
    stop(sprintf("Column '%s' (%s) has %d value(s) that are not numbers, e.g. %s (data rows %s). Fix them in a copy of the file or leave the cell empty.%s",
                 col, what, length(bad), paste(sprintf("'%s'", utils::head(unique(s[bad]), 4)), collapse = ", "),
                 paste(utils::head(bad, 6), collapse = ", "), hint), call. = FALSE)
  }
  v
}

# Bayesian point estimates: posterior medians (Meyer & Millar 1999 and common
# stock-assessment practice; also what the figures plot). Means of the
# right-skewed posteriors of K, steepness and the SDs are pulled up by the tail.
fs_use_medians <- function(est, chains) {
  all <- do.call(rbind, lapply(chains, as.matrix))
  med <- apply(all, 2, stats::median, na.rm = TRUE)
  est$estimate <- unname(med[est$param])
  est
}

# ---- figure helpers -------------------------------------------------------

# Entity colours from the user's palette (figures: palette / greyscale), in a
# fixed order so an entity (an estimator, the median trajectory) keeps its
# colour across every figure. For Okabe-Ito the high-contrast hues (blue,
# vermillion, green) come first; other palettes are used in their own order.
fs_pal <- function(n, opts) {
  p <- opts$palette %||% "okabe_ito"
  pal <- if (length(p) > 1) unlist(p) else eco_pal[[p]]
  if (identical(p, "okabe_ito")) pal <- pal[c(5, 6, 3, 1, 7, 2, 4, 8, 9)]
  pal[(seq_len(n) - 1) %% length(pal) + 1]
}

# Legend placement that survives the core theme (which replaces any theme()
# set here): figures with wide legends default to the bottom unless the user
# chose a position other than the default.
fs_legend_pos <- function(opts, default = "bottom") {
  l <- opts$legend %||% "right"
  if (identical(l, "right")) default else l
}

fs_inner <- function(opts) opts$inner_ci %||% 0.5

# Relabel parameter columns of a list of chain matrices (figures: labels:).
FS_PARAM_LABELS <- c(ln_alpha = "ln(alpha): productivity", beta = "beta: density dependence",
                     phi = "phi: residual autocorrelation", sigma = "sigma: recruitment SD (log)",
                     sigma_nu = "sigma_nu: innovation SD (log)", r = "r: intrinsic growth (per year)",
                     K = "K: carrying capacity", q = "q: catchability")
fs_relabel_chains <- function(chains, opts) lapply(chains, function(m) {
  m <- as.matrix(m)
  colnames(m) <- vapply(colnames(m), function(k) fig_label(k, if (k %in% names(FS_PARAM_LABELS)) FS_PARAM_LABELS[[k]] else k, opts), "")
  m })

# Trace / posterior panels: few x ticks, double width (tick labels of
# 4-5 small facets overlap at single-column width).
fs_mcmc_fig <- function(p, name, opts, aspect = 0.5)
  fs_size(p + ggplot2::scale_x_continuous(n.breaks = 4), name, opts, aspect = aspect, width = "double")

# Size hint in mm. The core uses a module's width/height hints but not an
# aspect hint, so convert aspect to a height here, unless the user set the
# figure's size in per_figure (then the user's values win).
fs_size <- function(p, name, opts, aspect, width = NULL) {
  u <- opts$per_figure[[name]]
  if (!is.null(u$height) || !is.null(u$aspect)) { attr(p, "fig_size") <- list(width = width); return(p) }
  w <- u$width %||% width %||% opts$width %||% "single"
  attr(p, "fig_size") <- list(width = width, height = fig_width_mm(w, opts) * aspect)
  p
}
