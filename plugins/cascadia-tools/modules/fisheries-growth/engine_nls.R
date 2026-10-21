# Likelihood engine: von Bertalanffy by nonlinear least squares (stats::nls,
# port algorithm), with starting values from FSA::findGrowthStarts (the
# current name of FSA::vbStarts), and Kimura's (1980) likelihood-ratio tests
# for differences in growth among groups.
#
# The parameters are estimated on the scale log Linf, log K, t0, so Wald
# intervals for Linf and K are asymmetric and stay positive. Error:
#   multiplicative: log(L) = log(VB(age)) + e, e ~ N(0, sigma)
#                   (the curve is the median length at age)
#   additive:       L = VB(age) + e,           e ~ N(0, sigma)  (Kimura 1980)

fg_starts_one <- function(d) {
  s <- tryCatch(suppressWarnings(suppressMessages(
         FSA::findGrowthStarts(length ~ age, data = d, type = "von Bertalanffy", param = 1))),
         error = function(e) NULL)
  ok <- !is.null(s) && all(is.finite(s)) && s[["Linf"]] > 0 && s[["K"]] > 0
  if (!ok) s <- c(Linf = 1.1 * max(d$length), K = 0.3, t0 = 0)
  s <- c(Linf = unname(s[["Linf"]]), K = unname(s[["K"]]), t0 = unname(s[["t0"]]))
  s[["K"]] <- min(max(s[["K"]], 0.01), 3)
  s[["Linf"]] <- max(s[["Linf"]], max(d$length) * 0.8)
  s[["t0"]] <- min(s[["t0"]], min(d$age) - 0.2)
  s
}

# Fit one VB model; `common` lists parameters shared by all groups.
fg_nls_one <- function(df, common = character(), error = "multiplicative", starts = NULL) {
  G <- nlevels(df$group); gi <- as.integer(df$group)
  if (is.null(starts)) starts <- t(vapply(split(df, df$group), fg_starts_one, numeric(3)))
  term <- function(p, lp) if (p %in% common || G == 1) lp else paste0(lp, "[gi]")
  pred <- sprintf("exp(%s) * (1 - exp(-exp(%s) * (age - %s)))",
                  term("Linf", "lLinf"), term("K", "lK"), term("t0", "t0"))
  f <- stats::as.formula(if (error == "multiplicative") paste("log(length) ~ log(", pred, ")")
                         else paste("length ~", pred))
  mk <- function(v, p) if (p %in% common || G == 1) mean(v) else unname(v)
  st <- list(lLinf = mk(log(starts[, "Linf"]), "Linf"), lK = mk(log(starts[, "K"]), "K"),
             t0 = mk(starts[, "t0"], "t0"))
  mina <- vapply(split(df$age, df$group), min, 0)
  t0up <- if ("t0" %in% common || G == 1) min(mina) - 0.01 else unname(mina) - 0.01
  st$t0 <- pmin(st$t0, t0up - 0.1)
  lo <- c(rep(log(1e-6), length(st$lLinf)), rep(log(1e-4), length(st$lK)), rep(-30, length(st$t0)))
  up <- c(rep(log(1e7), length(st$lLinf)), rep(log(20), length(st$lK)), t0up)
  dd <- data.frame(age = df$age, length = df$length, gi = gi)
  try_fit <- function(s) tryCatch(
    stats::nls(f, data = dd, start = s, algorithm = "port", lower = lo, upper = up,
               control = stats::nls.control(maxiter = 500, warnOnly = FALSE)),
    error = function(e) NULL)
  m <- try_fit(st)
  if (is.null(m)) {  # retry from perturbed starts
    for (j in 1:6) {
      s2 <- st; s2$lLinf <- s2$lLinf + c(0.3, -0.2, 0.6, 0, 0.15, -0.1)[j]
      s2$lK <- s2$lK + c(-0.5, 0.3, -1, 0.5, 0, -0.3)[j]
      s2$t0 <- pmin(s2$t0 - c(0.5, 0, 1, 0.2, 0, 0.5)[j], t0up - 0.1)
      m <- try_fit(s2); if (!is.null(m)) break
    }
  }
  m
}

# Kimura (1980) likelihood-ratio tests: general model (all three parameters
# group-specific) vs each constrained model. -2 ln(Lambda) = N ln(RSS_H/RSS_O).
fg_kimura <- function(df, error, general) {
  N <- nrow(df); G <- nlevels(df$group)
  rss_o <- sum(stats::residuals(general)^2)
  hyp <- list("Common Linf" = "Linf", "Common K" = "K", "Common t0" = "t0",
              "All parameters common" = c("Linf", "K", "t0"))
  rows <- lapply(names(hyp), function(h) {
    m <- fg_nls_one(df, hyp[[h]], error)
    if (is.null(m)) return(data.frame(hypothesis = h, rss = NA, chisq = NA, df = NA, p_value = NA))
    rss <- sum(stats::residuals(m)^2)
    x2 <- N * log(rss / rss_o); k <- length(hyp[[h]]) * (G - 1)
    data.frame(hypothesis = h, rss = rss, chisq = x2, df = k,
               p_value = stats::pchisq(x2, k, lower.tail = FALSE))
  })
  out <- do.call(rbind, rows)
  rbind(data.frame(hypothesis = "General (all group-specific)", rss = rss_o, chisq = NA, df = NA, p_value = NA), out)
}

# Gradient of the predicted (log-)length wrt (log Linf, log K, t0).
fg_grad <- function(age, lL, lK, t0, error) {
  Linf <- exp(lL); K <- exp(lK); e <- exp(-K * (age - t0)); pr <- Linf * (1 - e)
  if (error == "multiplicative")
    cbind(1, K * (age - t0) * e / (1 - e), -K * e / (1 - e))
  else cbind(pr, Linf * K * (age - t0) * e, -Linf * K * e)
}

# Parameter covariance per group from the general nls fit.
fg_group_par <- function(m, G) {
  co <- stats::coef(m); V <- stats::vcov(m)
  idx <- function(p, g) { nm <- if (G == 1) p else paste0(p, g); match(nm, names(co)) }
  lapply(seq_len(G), function(g) {
    ii <- c(idx("lLinf", g), idx("lK", g), idx("t0", g))
    list(theta = co[ii], V = V[ii, ii, drop = FALSE], ii = ii)
  })
}

fit_nls_fg <- function(data, settings) {
  df <- data$df; error <- settings$error
  G <- nlevels(df$group); gn <- levels(df$group); sg <- fg_safe(gn)
  z <- stats::qnorm(0.975)
  m <- fg_nls_one(df, character(), error)
  if (is.null(m)) stop("von Bertalanffy nls fit failed from every starting value: the data may not show curvature (too few old fish?)")
  co <- stats::coef(m); V <- stats::vcov(m)
  gp <- fg_group_par(m, G)
  rows <- list()
  for (g in seq_len(G)) {
    th <- gp[[g]]$theta; Vg <- gp[[g]]$V; se <- sqrt(pmax(diag(Vg), 0))
    K <- exp(th[2]); a90 <- th[3] + log(10) / K
    ga <- c(0, -log(10) / K, 1); se_a90 <- sqrt(drop(t(ga) %*% Vg %*% ga))
    rows[[length(rows) + 1]] <- data.frame(
      param = paste0(c("Linf_", "K_", "t0_", "age_at_90pct_Linf_"), sg[g]),
      estimate = c(exp(th[1]), K, th[3], a90),
      sd = c(exp(th[1]) * se[1], K * se[2], se[3], se_a90),
      lower = c(exp(th[1] - z * se[1]), exp(th[2] - z * se[2]), th[3] - z * se[3], a90 - z * se_a90),
      upper = c(exp(th[1] + z * se[1]), exp(th[2] + z * se[2]), th[3] + z * se[3], a90 + z * se_a90))
  }
  # Differences from the reference group (delta method on the natural scale)
  if (G > 1) {
    r <- 1
    for (g in setdiff(seq_len(G), r)) {
      for (p in 1:3) {
        gr <- numeric(length(co))
        nat <- function(k, gg) if (k < 3) exp(gp[[gg]]$theta[k]) else gp[[gg]]$theta[k]
        dv <- function(k, gg) if (k < 3) exp(gp[[gg]]$theta[k]) else 1
        gr[gp[[g]]$ii[p]] <- dv(p, g); gr[gp[[r]]$ii[p]] <- -dv(p, r)
        est <- nat(p, g) - nat(p, r); se <- sqrt(drop(t(gr) %*% V %*% gr))
        rows[[length(rows) + 1]] <- data.frame(
          param = paste0(c("Linf", "K", "t0")[p], "_diff_", sg[g], "_vs_", sg[r]),
          estimate = est, sd = se, lower = est - z * se, upper = est + z * se)
      }
    }
  }
  N <- nrow(df); rdf <- N - length(co); rss <- sum(stats::residuals(m)^2); s <- sqrt(rss / rdf)
  rows[[length(rows) + 1]] <- data.frame(param = "residual_sd", estimate = s, sd = s / sqrt(2 * rdf),
    lower = sqrt(rss / stats::qchisq(0.975, rdf)), upper = sqrt(rss / stats::qchisq(0.025, rdf)))
  est <- do.call(rbind, rows); rownames(est) <- NULL
  kim <- if (G > 1 && !isFALSE(settings$kimura)) fg_kimura(df, error, m) else NULL
  notes <- character()
  conv <- m$convInfo$isConv && all(is.finite(diag(V))) && all(diag(V) > 0)
  if (!conv) notes <- c(notes, "nls did not converge cleanly or the covariance matrix is singular.")
  list(estimates = est, diagnostics = list(converged = conv, notes = notes),
       nls = m, kimura = kim, gp = gp, sigma = s, rdf = rdf)
}

# Curve predictions (median length at age for multiplicative error) with a
# confidence band for the curve and a prediction interval for single fish.
fg_curve_nls <- function(fit, df, error, ci = 0.95, ages = NULL) {
  G <- nlevels(df$group); gn <- levels(df$group); z <- stats::qnorm(1 - (1 - ci) / 2)
  out <- lapply(seq_len(G), function(g) {
    a <- if (is.null(ages)) seq(max(0, min(df$age) - 0.5), max(df$age), length.out = 80) else ages
    th <- fit$gp[[g]]$theta; V <- fit$gp[[g]]$V
    t0 <- th[3]; a <- a[a > t0 + 0.05]
    J <- fg_grad(a, th[1], th[2], t0, error)
    se <- sqrt(pmax(rowSums((J %*% V) * J), 0))
    pr <- fg_vb(a, exp(th[1]), exp(th[2]), t0)
    sp <- sqrt(se^2 + fit$sigma^2)
    if (error == "multiplicative")
      data.frame(group = gn[g], x = a, estimate = pr, lower = pr * exp(-z * se), upper = pr * exp(z * se),
                 pi_lower = pr * exp(-z * sp), pi_upper = pr * exp(z * sp))
    else
      data.frame(group = gn[g], x = a, estimate = pr, lower = pr - z * se, upper = pr + z * se,
                 pi_lower = pr - z * sp, pi_upper = pr + z * sp)
  })
  do.call(rbind, out)
}
