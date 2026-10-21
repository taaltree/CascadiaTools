# Length-weight and condition.
#
# log W = log a + b log L + log(condition_g) + e,  e ~ N(0, s)
# fitted by least squares on the log-log scale with a common slope among
# groups (ANCOVA). The interaction model (group-specific b) is fitted too and
# compared by an F test, so a researcher knows whether a common b is
# defensible. Back-transformed predictions of MEAN weight carry the lognormal
# bias correction exp(s^2/2) (Sprugel 1983).
#
# Condition:
#   condition_ratio_<g>_vs_<ref> = exp(group effect): the weight of a fish of
#     group g relative to a reference-group fish of the same length.
#   Relative weight Wr = 100 W / Ws(L) when a standard-weight equation is
#   supplied (Wege & Anderson 1978; equations from FSA::wsVal / WSlit).
#   Otherwise relative condition Kn = W / W_hat(L) from the pooled fit
#   (Le Cren 1951).

fg_safe <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))

fit_lw <- function(data) {
  df <- data$df[is.finite(data$df$weight) & data$df$weight > 0, ]
  if (nrow(df) < 10) return(NULL)
  df$group <- droplevels(df$group)
  G <- nlevels(df$group); gn <- levels(df$group); sg <- fg_safe(gn)
  df$logL <- log(df$length); df$logW <- log(df$weight)
  m <- if (G > 1) stats::lm(logW ~ logL + group, data = df) else stats::lm(logW ~ logL, data = df)
  co <- stats::coef(m); V <- stats::vcov(m); rdf <- m$df.residual
  s <- summary(m)$sigma; tq <- stats::qt(0.975, rdf)
  se <- sqrt(diag(V))
  lref <- data$reference_length
  x <- c(1, log(lref), rep(0, length(co) - 2))
  lp <- sum(x * co); se_lp <- sqrt(drop(t(x) %*% V %*% x)); bc <- exp(s^2 / 2)
  rows <- data.frame(
    param = c("lw_b", "lw_a", "lw_residual_sd", paste0("mean_weight_at_", lref, "_", sg[1])),
    estimate = c(co[2], exp(co[1]), s, exp(lp) * bc),
    sd = c(se[2], exp(co[1]) * se[1], s / sqrt(2 * rdf), exp(lp) * bc * se_lp),
    lower = c(co[2] - tq * se[2], exp(co[1] - tq * se[1]), sqrt(s^2 * rdf / stats::qchisq(0.975, rdf)),
              exp(lp - tq * se_lp) * bc),
    upper = c(co[2] + tq * se[2], exp(co[1] + tq * se[1]), sqrt(s^2 * rdf / stats::qchisq(0.025, rdf)),
              exp(lp + tq * se_lp) * bc))
  slope_test <- NULL
  if (G > 1) {
    k <- paste0("group", gn[-1])
    rows <- rbind(rows, data.frame(param = paste0("condition_ratio_", sg[-1], "_vs_", sg[1]),
                                   estimate = exp(co[k]), sd = exp(co[k]) * se[k],
                                   lower = exp(co[k] - tq * se[k]), upper = exp(co[k] + tq * se[k])))
    mi <- stats::lm(logW ~ logL * group, data = df)
    av <- stats::anova(m, mi)
    slope_test <- list(F = av$F[2], df1 = av$Df[2], df2 = av$Res.Df[2], p_value = av$`Pr(>F)`[2],
                       b_by_group = vapply(gn, function(g) stats::coef(stats::lm(logW ~ logL, df[df$group == g, ]))[2], 0))
  }
  rownames(rows) <- NULL
  # Condition indices per fish
  ws <- data$ws
  if (!is.null(ws)) {
    df$Ws <- 10^(ws$log10_a + ws$b * log10(df$length))
    df$cond_index <- 100 * df$weight / df$Ws
    if (!is.null(ws$min_length)) df$cond_index[df$length < ws$min_length] <- NA
    idx_name <- "Relative weight Wr"
  } else {
    mp <- stats::lm(logW ~ logL, data = df)
    df$cond_index <- df$weight / exp(stats::fitted(mp) + summary(mp)$sigma^2 / 2)
    idx_name <- "Relative condition Kn"
  }
  cond_tab <- do.call(rbind, lapply(gn, function(g) {
    v <- df$cond_index[df$group == g]; v <- v[is.finite(v)]
    n <- length(v); mn <- mean(v); se <- stats::sd(v) / sqrt(n); tq <- if (n > 1) stats::qt(0.975, n - 1) else NA
    data.frame(group = g, n = n, mean = mn, lower = mn - tq * se, upper = mn + tq * se)
  }))
  list(lm = m, rows = rows, sigma = s, slope_test = slope_test, df = df,
       cond_tab = cond_tab, cond_name = idx_name, outliers = df$fish_id[abs(stats::rstandard(m)) > 4])
}
