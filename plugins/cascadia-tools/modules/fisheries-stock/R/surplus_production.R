# Bayesian state-space Schaefer surplus-production model (Meyer & Millar 1999).
#
#   P_t = B_t / K
#   P_1 ~ lognormal(log(p1_median), p1_sd)
#   P_t ~ lognormal(log(max(P_{t-1} + r P_{t-1}(1 - P_{t-1}) - C_{t-1}/K, 0.001)), sigma_proc)
#   I_t ~ lognormal(log(q K P_t), sigma_obs)
#   MSY = rK/4, B_MSY = K/2, F_MSY (exploitation rate) = r/2
#   B/B_MSY = 2 P_t, F/F_MSY = (C_t / B_t) / (r/2)
#
# Priors (defaults, all overridable in settings):
#   r ~ lognormal(log(prior_r[1]), prior_r[2])        default median 0.3, sdlog 0.6
#   K ~ lognormal(log(prior_K[1]), prior_K[2])        default median 10 x mean catch, sdlog 1
#   log q ~ uniform (wide; ~ MM99's 1/q gamma(0.001, 0.001))
#   1/sigma_proc^2 ~ gamma(3.785518, 0.010223)        Meyer & Millar (1999)
#   1/sigma_obs^2  ~ gamma(1.708603, 0.008613854)     Meyer & Millar (1999)

sp_simulate <- function(params, design, seed) {
  set.seed(seed)
  Tn <- design$sp_years %||% 30
  r <- params$r; K <- params$K; q <- params$q
  sp <- params$sigma_proc; so <- params$sigma_obs
  pattern <- design$catch_pattern %||% "contrast"
  # Exploitation relative to U_MSY = r/2 over time
  t01 <- (seq_len(Tn) - 1) / (Tn - 1)
  fpk <- design$F_peak %||% 2
  f_rel <- switch(pattern,
    contrast = ifelse(t01 <= 0.45, 0.2 + (fpk - 0.2) * t01 / 0.45,
                      ifelse(t01 <= 0.65, fpk + (0.5 - fpk) * (t01 - 0.45) / 0.2, 0.5)),
    one_way = 0.2 + ((design$F_end %||% 1.6) - 0.2) * t01,
    stop("catch_pattern must be 'contrast' or 'one_way'"))
  f_rel <- f_rel * exp(stats::rnorm(Tn, 0, design$F_cv %||% 0.1))
  u <- pmin(f_rel * r / 2, 0.9)
  B <- numeric(Tn); C <- numeric(Tn)
  B[1] <- K * (params$P1 %||% 1) * exp(stats::rnorm(1, 0, sp))
  for (t in seq_len(Tn)) {
    C[t] <- u[t] * B[t]
    if (t < Tn) B[t + 1] <- max(B[t] + r * B[t] * (1 - B[t] / K) - C[t], 0.001 * K) * exp(stats::rnorm(1, 0, sp))
  }
  I <- q * B * exp(stats::rnorm(Tn, 0, so))
  miss <- design$index_missing %||% 0
  if (miss > 0) I[sample(seq_len(Tn), round(miss * Tn))] <- NA
  df <- data.frame(year = 1990 + seq_len(Tn), catch = C, index = I)
  truth <- c(intrinsic_growth_r = r, carrying_capacity_K = K, catchability_q = q, MSY = r * K / 4,
             B_over_BMSY_final = 2 * B[Tn] / K, F_over_FMSY_final = (C[Tn] / B[Tn]) / (r / 2),
             depletion_final = B[Tn] / K)
  if (isTRUE(design$truth_sigmas)) truth <- c(truth, sigma_process = sp, sigma_observation = so)
  list(data = list(model = "surplus_production", df = df,
                   meta = list(source = "simulated", catch_units = "tonnes", index_name = "CPUE", B_true = B)),
       truth = truth)
}

sp_code <- function() nimble::nimbleCode({
  logr ~ dnorm(mu_r, sd = sd_r)
  logK ~ dnorm(mu_K, sd = sd_K)
  logq ~ dunif(lq_lo, lq_hi)
  isigma2 ~ dgamma(a0, b0)
  itau2 ~ dgamma(c0, d0)
  r <- exp(logr); K <- exp(logK); q <- exp(logq)
  P[1] ~ dlnorm(log(p1_median), sdlog = p1_sd)
  for (t in 2:nT) {
    Pm[t] <- max(P[t - 1] + r * P[t - 1] * (1 - P[t - 1]) - C[t - 1] / K, 0.001)
    P[t] ~ dlnorm(log(Pm[t]), tau = isigma2)
  }
  for (i in 1:nI) {
    I[i] ~ dlnorm(logq + logK + log(P[iyr[i]]), tau = itau2)
  }
})

sp_mcmc <- function(settings) {
  fast <- isTRUE(settings$fast)
  list(niter = settings$niter %||% if (fast) 30000 else 100000, nburnin = settings$nburnin %||% if (fast) 10000 else 30000,
       thin = settings$thin %||% if (fast) 10 else 20)
}

sp_fit <- function(data, settings) {
  suppressPackageStartupMessages(library(nimble))  # model building needs nimble attached
  df <- data$df
  Tn <- nrow(df); C <- df$catch
  obs <- which(!is.na(df$index) & df$index > 0)
  I <- df$index[obs]
  num <- function(x) if (is.null(x)) NULL else as.numeric(unlist(x))
  pr <- num(settings$prior_r) %||% c(0.3, 0.6)
  pK <- num(settings$prior_K) %||% c(10 * mean(C), 1)
  ig <- num(settings$prior_sigma_proc) %||% c(3.785518, 0.010223)
  io <- num(settings$prior_sigma_obs) %||% c(1.708603, 0.008613854)
  p1 <- num(settings$prior_P1) %||% c(1, 0.2)
  lq_mid <- log(mean(I)) - log(pK[1])
  consts <- list(nT = Tn, nI = length(obs), iyr = obs)
  dat <- list(C = C, I = I, mu_r = log(pr[1]), sd_r = pr[2], mu_K = log(pK[1]), sd_K = pK[2],
              lq_lo = lq_mid - 15, lq_hi = lq_mid + 15, a0 = ig[1], b0 = ig[2], c0 = io[1], d0 = io[2],
              p1_median = p1[1], p1_sd = p1[2])
  # Inits: plausible r, K from the prior; P from the scaled index
  inits <- function() {
    r0 <- pr[1] * exp(stats::rnorm(1, 0, 0.2)); K0 <- max(pK[1] * exp(stats::rnorm(1, 0, 0.2)), 2 * max(C))
    P0 <- numeric(Tn); P0[1] <- p1[1]
    for (t in 2:Tn) P0[t] <- max(P0[t - 1] + r0 * P0[t - 1] * (1 - P0[t - 1]) - C[t - 1] / K0, 0.05)
    lq0 <- mean(log(I) - log(K0 * P0[obs]))
    list(logr = log(r0), logK = log(K0), logq = lq0, isigma2 = 100, itau2 = 50, P = P0)
  }
  mc <- sp_mcmc(settings)
  set.seed(settings$seed %||% 1)
  res <- nimble_run(sp_code(), constants = consts, data = dat, inits = inits,
                    monitors = c("r", "K", "q", "isigma2", "itau2", "P"),
                    niter = mc$niter, nburnin = mc$nburnin,
                    nchains = settings$nchains %||% 3, thin = mc$thin,
                    cache_key = paste("fs_sp", Tn, length(obs), paste(obs, collapse = "."), sep = "_"),
                    configure = function(conf, m) {
                      conf$removeSamplers(c("logr", "logK", "logq"))
                      conf$addSampler(c("logr", "logK", "logq"), "AF_slice")
                      conf$addSampler(c("logr", "logK"), "RW_block", control = list(adaptInterval = 100))
                    })
  pcols <- paste0("P[", seq_len(Tn), "]")
  chains <- lapply(res$samples, function(m) {
    r <- m[, "r"]; K <- m[, "K"]; P <- m[, pcols, drop = FALSE]
    PT <- P[, Tn]; BT <- K * PT
    cbind(intrinsic_growth_r = r, carrying_capacity_K = K, catchability_q = m[, "q"], MSY = r * K / 4,
          B_MSY = K / 2, F_MSY = r / 2,
          B_over_BMSY_final = 2 * PT, F_over_FMSY_final = (C[Tn] / BT) / (r / 2), depletion_final = PT,
          sigma_process = 1 / sqrt(m[, "isigma2"]), sigma_observation = 1 / sqrt(m[, "itau2"]))
  })
  est <- fs_use_medians(fs_summarise_chains(chains), chains)
  allm <- stack_draws(res$samples)
  P <- allm[, pcols, drop = FALSE]
  rh_all <- max(c(est$rhat, res$estimates$rhat[res$estimates$param %in% pcols]), na.rm = TRUE)
  es <- min(est$ess_bulk, na.rm = TRUE)
  traj <- list(year = df$year, P = P, r = allm[, "r"], K = allm[, "K"], q = allm[, "q"])
  notes <- sprintf("Bayesian state-space Schaefer model in Nimble: %d chains x %d iterations (%d burn-in, thin %d). Max R-hat %.3f (parameters and biomass states), min bulk ESS %.0f.",
                   settings$nchains %||% 3, mc$niter, mc$nburnin, mc$thin, rh_all, es)
  # prior-vs-posterior overlap for r and K: data informativeness
  shrink <- c(r = stats::sd(log(allm[, "r"])) / pr[2], K = stats::sd(log(allm[, "K"])) / pK[2])
  if (any(shrink > 0.7))
    notes <- c(notes, sprintf("The data add little information about %s (posterior SD on the log scale is %s of the prior SD): estimates are driven by the prior. Typical cause: a 'one-way trip' catch history without a period of stock recovery or low catch.",
                              paste(names(shrink)[shrink > 0.7], collapse = " and "), paste(sprintf("%.0f%%", 100 * shrink[shrink > 0.7]), collapse = " and ")))
  tr_cor <- stats::cor(log(allm[, "r"]), log(allm[, "K"]))
  list(estimates = est,
       diagnostics = list(converged = is.finite(rh_all) && rh_all < 1.05, rhat_max = rh_all, ess_min = es, notes = notes,
                          prior_posterior_sd_ratio = shrink, cor_logr_logK = tr_cor, t_total = res$diagnostics$t_total),
       object = list(model = "surplus_production", df = df, traj = traj, samples = res$samples, obs = obs,
                     priors = list(r = pr, K = pK, sigma_proc = ig, sigma_obs = io, P1 = p1)),
       settings = settings)
}

sp_read_data <- function(config, dc) {
  raw <- utils::read.csv(file.path(config$.dir, dc$file), check.names = FALSE, stringsAsFactors = FALSE)
  yc <- dc$year_col %||% "year"; cc <- dc$catch_col %||% "catch"; ic <- dc$index_col %||% "cpue"
  df <- data.frame(year = fs_num(fs_col(raw, yc, "year"), yc, "year"),
                   catch = fs_num(fs_col(raw, cc, "catch"), cc, "catch"),
                   index = fs_num(fs_col(raw, ic, "index"), ic, "index"))
  df <- df[order(df$year), ]
  list(model = "surplus_production", df = df,
       meta = list(file = dc$file, species = dc$species %||% "stock", catch_units = dc$catch_units %||% "tonnes",
                   index_name = dc$index_name %||% "CPUE"))
}

sp_check_data <- function(data) {
  pr <- character(); df <- data$df
  if (anyNA(df$catch)) pr <- c(pr, sprintf("Catch is missing in years %s: the model needs a catch for every year (use 0 or an estimate).", paste(df$year[is.na(df$catch)], collapse = ", ")))
  if (any(df$catch < 0, na.rm = TRUE)) pr <- c(pr, "Negative catches.")
  if (any(diff(df$year) != 1)) {
    gap <- setdiff(seq(min(df$year), max(df$year)), df$year); dup <- unique(df$year[duplicated(df$year)])
    pr <- c(pr, sprintf("Years must be consecutive, one row per year%s%s; add rows with the catch and an empty index for gap years.",
                        if (length(gap)) paste0(" (missing: ", paste(gap, collapse = ", "), ")") else "",
                        if (length(dup)) paste0(" (repeated: ", paste(dup, collapse = ", "), ")") else ""))
  }
  if (any(df$index <= 0, na.rm = TRUE)) pr <- c(pr, sprintf("Index values must be > 0 (years %s); leave the cell empty for years without a survey.", paste(df$year[!is.na(df$index) & df$index <= 0], collapse = ", ")))
  if (sum(!is.na(df$index)) < 8) pr <- c(pr, sprintf("Only %d years with an abundance index: a production model needs ~10 or more.", sum(!is.na(df$index))))
  pr
}

sp_kobe_status <- function(b, f) {
  ifelse(b >= 1 & f <= 1, "not overfished, no overfishing",
         ifelse(b >= 1 & f > 1, "not overfished, but overfishing",
                ifelse(b < 1 & f <= 1, "overfished, no overfishing", "overfished and overfishing")))
}

sp_interpret <- function(fit, data) {
  e <- fit$estimates; g <- function(p) e[e$param == p, ]
  o <- fit$object; cu <- data$meta$catch_units %||% "tonnes"
  Tn <- nrow(o$df); yr <- o$df$year[Tn]
  ch <- stack_draws(lapply(o$samples, as.matrix))
  bb <- 2 * ch[, paste0("P[", Tn, "]")]; ff <- (o$df$catch[Tn] / (ch[, "K"] * ch[, paste0("P[", Tn, "]")])) / (ch[, "r"] / 2)
  st <- table(factor(sp_kobe_status(bb, ff), levels = c("not overfished, no overfishing", "not overfished, but overfishing", "overfished, no overfishing", "overfished and overfishing"))) / length(bb)
  r <- g("intrinsic_growth_r"); K <- g("carrying_capacity_K"); M <- g("MSY"); B <- g("B_over_BMSY_final"); F <- g("F_over_FMSY_final")
  out <- c(
    sprintf("The stock's maximum sustainable yield is estimated at %s %s per year (95%% credible interval %s to %s), from an intrinsic growth rate r = %.2f (%.2f to %.2f) and carrying capacity K = %s %s (%s to %s). The average catch over the last 5 years was %s %s.",
            fs_sig(M$estimate), cu, fs_sig(M$lower), fs_sig(M$upper), r$estimate, r$lower, r$upper, fs_sig(K$estimate), cu, fs_sig(K$lower), fs_sig(K$upper),
            fs_sig(mean(utils::tail(o$df$catch, 5))), cu),
    sprintf("In %s, biomass was %.2f times the biomass that supports MSY (B/B_MSY, %.2f to %.2f) and the harvest rate was %.2f times the MSY rate (F/F_MSY, %.2f to %.2f). The posterior probability that the stock is overfished (B < B_MSY) is %s and that overfishing is occurring (F > F_MSY) is %s.",
            yr, B$estimate, B$lower, B$upper, F$estimate, F$lower, F$upper, fs_pct(mean(bb < 1)), fs_pct(mean(ff > 1))),
    sprintf("Kobe-plot status probabilities for %s: %s.", yr, paste(sprintf("%s %s", names(st), fs_pct(as.numeric(st))), collapse = "; ")))
  sh <- fit$diagnostics$prior_posterior_sd_ratio
  if (!is.null(sh) && any(sh > 0.7))
    out <- c(out, "Caution: the catch and index series carry little information about productivity and stock size (posterior close to prior). Status estimates then reflect the priors more than the data; report them as conditional on the priors and seek contrast (a period of reduced catch and recovery) or informative priors from related stocks.")
  c(out, sprintf("r and K are strongly negatively correlated in the posterior (correlation of logs %.2f): many combinations fit the data, which is why MSY and relative status (B/B_MSY) are usually better determined than r or K alone. The model assumes the index is proportional to biomass (constant catchability) and that catches are known without error.", fit$diagnostics$cor_logr_logK))
}

sp_methods <- function(fit, data, settings) {
  o <- fit$object; p <- o$priors
  sprintf(paste(
    "Stock dynamics were modelled with a Bayesian state-space Schaefer surplus-production model (Schaefer 1954; Meyer & Millar 1999) fitted to %d years of catch and %d years of an abundance index.",
    "Biomass relative to carrying capacity, P_t = B_t/K, followed P_t = (P_{t-1} + r P_{t-1}(1 - P_{t-1}) - C_{t-1}/K) exp(u_t) with lognormal process error (SD sigma_proc), and the index followed I_t = q K P_t exp(v_t) with lognormal observation error (SD sigma_obs).",
    "Priors: r ~ lognormal(median %.3g, log-SD %.2f); K ~ lognormal(median %.3g, log-SD %.2f); log q ~ uniform over a 30-log-unit range; 1/sigma_proc^2 ~ Gamma(%.4g, %.4g) and 1/sigma_obs^2 ~ Gamma(%.4g, %.4g) (Meyer & Millar 1999); P_1 ~ lognormal(median %.2f, log-SD %.2f).",
    "Reference points were MSY = rK/4, B_MSY = K/2 and F_MSY = r/2 (exploitation rate); stock status is B/B_MSY and F/F_MSY in the final year, shown in a Kobe plot.",
    "The model was fitted in Nimble (de Valpine et al. 2017) with %d chains of %d iterations (%d burn-in, thinning %d), with a joint automated-factor slice sampler on (log r, log K, log q); convergence was assessed with R-hat and bulk ESS. Point estimates are posterior medians with 95%% equal-tailed credible intervals."),
    nrow(o$df), length(o$obs), p$r[1], p$r[2], p$K[1], p$K[2], p$sigma_proc[1], p$sigma_proc[2], p$sigma_obs[1], p$sigma_obs[2], p$P1[1], p$P1[2],
    settings$nchains %||% 3, sp_mcmc(settings)$niter, sp_mcmc(settings)$nburnin, sp_mcmc(settings)$thin)
}

sp_figures <- function(fit, data, opts) {
  o <- fit$object; df <- o$df; tr <- o$traj; Tn <- nrow(df)
  ci <- opts$ci %||% 0.95; ici <- fs_inner(opts); a <- (1 - ci) / 2; ia <- (1 - ici) / 2
  cu <- data$meta$catch_units %||% "tonnes"
  col <- fs_pal(2, opts)                      # 1 = model estimate (median trajectory), 2 = reference point
  qf <- function(M, p) apply(M, 2, stats::quantile, p, na.rm = TRUE)
  bands <- function(M) data.frame(estimate = qf(M, 0.5), lower = qf(M, a), upper = qf(M, 1 - a), lower50 = qf(M, ia), upper50 = qf(M, 1 - ia))
  band_layers <- function(colr) list(
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = colr, alpha = 0.15),
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = colr, alpha = 0.25),
    ggplot2::geom_line(colour = colr, linewidth = 0.8))
  lab_year <- fig_label("year", "Year", opts)
  lab_B <- fig_label("biomass", sprintf("Biomass (%s)", cu), opts)
  B <- tr$P * tr$K
  bt <- data.frame(time = df$year, bands(B))
  qm <- stats::median(tr$q); bmsy <- stats::median(tr$K) / 2
  obs <- data.frame(time = df$year[o$obs], y = df$index[o$obs] / qm)
  pb <- ggplot2::ggplot(bt, ggplot2::aes(time, estimate)) +
    ggplot2::geom_col(data = df, ggplot2::aes(year, catch), inherit.aes = FALSE, fill = "#bdbdbd", width = 0.7) +
    band_layers(col[1]) +
    ggplot2::geom_hline(yintercept = bmsy, linetype = 2, colour = col[2], linewidth = 0.5) +
    ggplot2::annotate("text", x = max(df$year), y = bmsy, label = fig_label("B_MSY_lab", "B[MSY]", opts), parse = TRUE,
                      hjust = 1, vjust = -0.4, size = 2.5, colour = col[2]) +
    ggplot2::geom_point(data = obs, ggplot2::aes(time, y), inherit.aes = FALSE, size = 1.1, shape = 21, fill = "white", colour = eco_pal$ink, stroke = 0.5) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.04)), labels = scales::label_comma()) +
    ggplot2::labs(x = lab_year, y = fig_label("biomass_catch", sprintf("Biomass or catch (%s)", cu), opts))
  pb <- fs_size(pb, "biomass_trajectory", opts, aspect = 0.65)
  # Kobe: trajectory of medians + final-year posterior cloud
  Bm <- 2 * tr$P
  Fm <- sweep(1 / (tr$P * tr$K), 2, df$catch, "*") / (tr$r / 2)
  kt <- data.frame(year = df$year, b = qf(Bm, 0.5), f = qf(Fm, 0.5))
  set.seed(1)
  idx <- sample(nrow(Bm), min(1500, nrow(Bm)))
  cloud <- data.frame(b = Bm[idx, Tn], f = Fm[idx, Tn])
  xmax <- max(2, stats::quantile(Bm[, Tn], 0.99), max(kt$b)) * 1.05
  ymax <- max(2, stats::quantile(Fm[, Tn], 0.99), max(kt$f)) * 1.05
  qlab_txt <- c(fig_label("kobe_red", "Overfished,\noverfishing", opts), fig_label("kobe_green", "Not overfished,\nno overfishing", opts),
                fig_label("kobe_yellow", "Overfished,\nno overfishing", opts), fig_label("kobe_orange", "Not overfished,\noverfishing", opts))
  # Kobe colours are a management convention (red/yellow/orange/green); greys when printing in greyscale.
  qfill <- if (isTRUE(opts$greyscale) || identical(opts$palette, "grey")) c("#8c8c8c", "#f2f2f2", "#c8c8c8", "#c8c8c8")
           else c("#D55E00", "#009E73", "#F0E442", "#E69F00")
  quad <- data.frame(xmin = c(0, 1, 0, 1), xmax = c(1, Inf, 1, Inf), ymin = c(1, 0, 0, 1), ymax = c(Inf, 1, 1, Inf), fill = qfill)
  qlab <- data.frame(x = c(0.03, 0.97, 0.03, 0.97) * xmax, y = c(0.97, 0.03, 0.03, 0.97) * ymax,
                     hj = c(0, 1, 0, 1), vj = c(1, 0, 0, 1), lab = qlab_txt)
  ends <- kt[c(1, Tn), ]
  pk <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = quad, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = quad$fill, alpha = 0.3) +
    ggplot2::geom_text(data = qlab, ggplot2::aes(x, y, label = lab, hjust = hj, vjust = vj), size = 2.2, colour = eco_pal$ink, lineheight = 0.9) +
    ggplot2::geom_point(data = cloud, ggplot2::aes(b, f), size = 0.35, alpha = 0.25, colour = eco_pal$muted) +
    ggplot2::geom_hline(yintercept = 1, linewidth = 0.3) + ggplot2::geom_vline(xintercept = 1, linewidth = 0.3) +
    ggplot2::geom_path(data = kt, ggplot2::aes(b, f), colour = col[1], linewidth = 0.5) +
    ggplot2::geom_point(data = kt, ggplot2::aes(b, f), colour = col[1], size = 0.8) +
    ggplot2::geom_point(data = ends, ggplot2::aes(b, f), shape = 21, fill = "white", colour = col[1], size = 2.2, stroke = 0.8)
  pk <- pk + if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_label_repel(data = ends, ggplot2::aes(b, f, label = year), size = 2.3, colour = eco_pal$ink, fill = "white",
                              label.size = 0, label.padding = 0.1, box.padding = 0.3, min.segment.length = 0, segment.size = 0.25, seed = 1,
                              # push year labels toward the panel centre, away from the corner quadrant labels
                              nudge_x = sign(xmax / 2 - ends$b) * 0.09 * xmax, nudge_y = sign(ymax / 2 - ends$f) * 0.09 * ymax)
    else ggplot2::geom_text(data = ends, ggplot2::aes(b, f, label = year), size = 2.3, vjust = -1.1, colour = eco_pal$ink)
  pk <- pk + ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    ggplot2::labs(x = fig_label("B_over_BMSY", "Relative biomass, B / B_MSY", opts), y = fig_label("F_over_FMSY", "Relative fishing mortality, F / F_MSY", opts))
  pk <- fs_size(pk, "kobe", opts, aspect = 1)
  # Production curve with bands + "observed" surplus production from the median biomass
  Bg <- seq(0, max(stats::median(tr$K), max(bt$estimate)) * 1.05, length.out = 120)
  ns <- min(2000, length(tr$r)); ii <- round(seq(1, length(tr$r), length.out = ns))
  SPm <- vapply(Bg, function(b) tr$r[ii] * b * (1 - b / tr$K[ii]), numeric(ns))
  spd <- data.frame(x = Bg, bands(SPm))
  bmed <- bt$estimate
  osp <- data.frame(x = bmed[-Tn], y = diff(bmed) + df$catch[-Tn], year = df$year[-Tn])
  pp <- ggplot2::ggplot(spd, ggplot2::aes(x, estimate)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = eco_pal$null) +
    band_layers(col[1]) +
    ggplot2::geom_path(data = osp, ggplot2::aes(x, y), inherit.aes = FALSE, colour = "#b3b3b3", linewidth = 0.3) +
    ggplot2::geom_point(data = osp, ggplot2::aes(x, y, fill = year), inherit.aes = FALSE, size = 1.4, shape = 21, colour = eco_pal$ink, stroke = 0.2) +
    ggplot2::scale_fill_viridis_c(name = lab_year, option = "viridis",
                                  guide = ggplot2::guide_colourbar(barheight = ggplot2::unit(0.4, "lines"), barwidth = ggplot2::unit(8, "lines"),
                                                                   position = fs_legend_pos(opts), title.vjust = 0.9)) +
    ggplot2::scale_x_continuous(labels = scales::label_comma()) + ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::coord_cartesian(ylim = c(min(0, osp$y), max(spd$upper, osp$y))) +
    ggplot2::labs(x = lab_B, y = fig_label("surplus", sprintf("Surplus production (%s/year)", cu), opts))
  pp <- fs_size(pp, "surplus_production", opts, aspect = 0.85)
  # Index fit
  Ipred <- sweep(B, 1, tr$q, "*")
  it <- data.frame(time = df$year, bands(Ipred))
  pi <- ggplot2::ggplot(it, ggplot2::aes(time, estimate)) + band_layers(col[1]) +
    ggplot2::geom_point(data = data.frame(time = df$year[o$obs], y = df$index[o$obs]), ggplot2::aes(time, y), inherit.aes = FALSE,
                        size = 1.2, colour = eco_pal$ink) +
    ggplot2::labs(x = lab_year, y = fig_label("index", data$meta$index_name %||% "Abundance index", opts))
  pi <- fs_size(pi, "index_fit", opts, aspect = 0.6)
  ch <- lapply(o$samples, function(m) m[, c("r", "K", "q"), drop = FALSE])
  dr <- stack_draws(lapply(o$samples, as.matrix))
  pl <- c(r = fig_label("log10_r", "log10 r (per year)", opts), K = fig_label("log10_K", sprintf("log10 K (%s)", cu), opts),
          MSY = fig_label("log10_MSY", sprintf("log10 MSY (%s/year)", cu), opts))
  bnd <- sprintf("%d%% and %d%% credible intervals", round(100 * ci), round(100 * ici))
  list(
    biomass_trajectory = with_caption(pb, sprintf(paste(
      "Estimated biomass (line: posterior median; bands: %s); open points: abundance index divided by the median catchability q;",
      "grey bars: catch; dashed line: B_MSY = K/2 (posterior median, %s %s)."), bnd, fs_sig(bmsy), cu)),
    kobe = with_caption(pk, sprintf(paste(
      "Kobe plot: posterior median trajectory of B/B_MSY and F/F_MSY (line; open circles mark %d and %d) and %d posterior draws for %d (grey cloud,",
      "the uncertainty in final-year status). Quadrants: green = not overfished and no overfishing; yellow = overfished, no overfishing;",
      "orange = overfishing, not overfished; red = overfished and overfishing."),
      df$year[1], df$year[Tn], nrow(cloud), df$year[Tn])),
    surplus_production = with_caption(pp, sprintf(paste(
      "Schaefer surplus-production curve rB(1 - B/K) (line: posterior median; bands: %s) with the annual surplus production implied",
      "by the median biomass trajectory (B_{t+1} - B_t + C_t; points shaded by year, joined in time order)."), bnd)),
    index_fit = with_caption(pi, sprintf("Observed abundance index (points) and the model-predicted index qB_t (line: posterior median; bands: %s).", bnd)),
    trace = with_caption(fs_mcmc_fig(fig_trace(fs_relabel_chains(ch, opts)), "trace", opts, aspect = 0.35), "MCMC traces for r, K and q (one colour per chain); chains should overlap with no trend."),
    posterior = with_caption({
      pr_r <- stats::rlnorm(20000, log(o$priors$r[1]), o$priors$r[2]); pr_K <- stats::rlnorm(20000, log(o$priors$K[1]), o$priors$K[2])
      post <- cbind(log10(dr[, "r"]), log10(dr[, "K"]), log10(dr[, "r"] * dr[, "K"] / 4)); colnames(post) <- pl
      prior <- cbind(log10(pr_r), log10(pr_K), log10(pr_r * pr_K / 4)); colnames(prior) <- pl
      pp2 <- fig_posterior(post, prior = prior, ncol = 3)
      fs_mcmc_fig(pp2, "posterior", opts, aspect = 0.35) },
      "Posterior densities (shaded) and prior densities (dashed) of r, K and MSY on a log10 scale. A posterior that matches its prior means the data are uninformative about that quantity.")
  )
}
