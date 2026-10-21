# jolly-seber: simulate open-population capture-recapture data with known truth.
#
# The generating process is the Schwarz & Arnason (1996) superpopulation
# model: Nsuper animals use the study area at some time during the study;
# each enters at occasion e with probability b_e (entry_probs, sum 1), is
# alive and present at e, and survives each interval with probability
# phi^dt (per time unit, optional group effect on the logit scale). Alive
# animals are captured with probability p_t. Optional Markovian temporary
# emigration (design$temp_emigration = c(leave, return)) makes alive animals
# unavailable for capture; it is NOT in the fitted model and serves as the
# negative control.

js_entry_probs <- function(params, T) {
  if (!is.null(params$entry_probs)) {
    b <- unlist(params$entry_probs)
    if (length(b) != T) stop("params$entry_probs must have one value per occasion (", T, ")")
    return(b / sum(b))
  }
  b1 <- params$entry_prob_first
  c(b1, rep((1 - b1) / (T - 1), T - 1))
}

js_simulate <- function(params, design, seed) {
  set.seed(seed)
  T <- as.integer(design$n_occasions)
  if (T < 3) stop("Jolly-Seber needs at least 3 occasions")
  dt <- rep(unlist(design$interval), length.out = T - 1)
  Ns <- as.integer(round(params$superpopulation_size))
  b <- js_entry_probs(params, T)
  # groups (e.g. sex)
  has_g <- isTRUE(design$groups)
  lev <- unlist(design$group_levels)
  g <- if (has_g) stats::rbinom(Ns, 1, params[[paste0("prop_", lev[2])]]) else rep(0L, Ns)
  bg <- if (has_g) params[[paste0("survival_beta_", lev[2])]] else 0
  # survival per time unit (constant or by interval) on the logit scale
  lphi <- if (!is.null(params$survival_by_interval)) stats::qlogis(rep(unlist(params$survival_by_interval), length.out = T - 1))
          else rep(stats::qlogis(params$survival), T - 1)
  # detection: constant, by occasion, or logit-linear in a standardised effort covariate
  has_cov <- isTRUE(design$detection_covariate)
  zc <- if (has_cov) as.numeric(scale(stats::runif(T))) else rep(0, T)
  lp <- if (!is.null(params$detection_by_occasion)) stats::qlogis(rep(unlist(params$detection_by_occasion), length.out = T))
        else stats::qlogis(params$detection) + (if (has_cov) params$detection_beta_effort * zc else 0)
  p <- stats::plogis(lp)
  # entry, survival
  e <- sample.int(T, Ns, replace = TRUE, prob = b)
  alive <- matrix(0L, Ns, T)
  for (t in seq_len(T)) {
    if (t == 1) alive[, 1] <- as.integer(e == 1)
    else {
      ph <- stats::plogis(lphi[t - 1] + bg * g)^dt[t - 1]
      surv <- alive[, t - 1] == 1 & stats::runif(Ns) < ph
      alive[, t] <- as.integer(surv | e == t)
    }
  }
  # availability (temporary emigration; negative control only)
  avail <- matrix(1L, Ns, T)
  te <- unlist(design$temp_emigration)
  if (length(te) == 2 && any(te > 0)) {
    for (t in seq_len(T)) {
      if (t == 1) next
      prev <- avail[, t - 1]
      stay_away <- stats::runif(Ns) >= te[2]   # away -> still away
      leave <- stats::runif(Ns) < te[1]        # present -> away
      avail[, t] <- ifelse(e == t, 1L, ifelse(prev == 1, as.integer(!leave), as.integer(!stay_away)))
    }
  }
  y <- matrix(stats::rbinom(Ns * T, 1, rep(p, each = Ns)), Ns, T) * alive * avail
  det <- rowSums(y) > 0
  ch <- y[det, , drop = FALSE]
  rownames(ch) <- sprintf("A%04d", seq_len(sum(det)))
  grp <- NULL
  if (has_g) {
    grp <- lev[g[det] + 1]
    pu <- if (is.null(design$prop_group_unknown)) 0 else design$prop_group_unknown
    if (pu > 0) grp[stats::runif(length(grp)) < pu] <- NA
  }
  yr0 <- if (is.null(design$first_year)) 1 else design$first_year
  occ <- data.frame(occasion = seq_len(T), label = as.character(yr0 + c(0, cumsum(dt))),
                    interval = c(dt, NA))
  if (has_cov) occ$effort <- round(design$effort_mean + design$effort_sd * zc, 1)
  data <- list(ch = ch, group = grp, group_levels = if (has_g) lev else NULL,
               group_name = if (has_g) (if (is.null(design$group_name)) "group" else design$group_name) else NULL,
               occasions = occ, detection_covariates = if (has_cov) "effort" else character(),
               time_unit = if (is.null(design$time_unit)) "year" else design$time_unit)
  # realised truth
  N <- colSums(alive)
  B <- vapply(2:T, function(t) sum(e == t), 1)
  tr <- c(superpopulation_size = Ns,
          stats::setNames(N, paste0("abundance_", seq_len(T))),
          stats::setNames(B, paste0("entrants_", 2:T)),
          stats::setNames(b, paste0("entry_prob_", seq_len(T))),
          stats::setNames(B / N[-T], paste0("recruitment_rate_", seq_len(T - 1))),
          stats::setNames(N[-1] / N[-T], paste0("growth_rate_", seq_len(T - 1))),
          growth_rate_mean = (N[T] / N[1])^(1 / sum(dt)))
  if (is.null(params$survival_by_interval)) tr <- c(tr, survival = params$survival)
  else tr <- c(tr, stats::setNames(stats::plogis(lphi), paste0("survival_", seq_len(T - 1))))
  if (has_g) {
    tr <- c(tr, stats::setNames(bg, paste0("survival_beta_", lev[2])),
            stats::setNames(mean(g), paste0("prop_", lev[2])))   # realised share in the superpopulation
    if (is.null(params$survival_by_interval))
      tr <- c(tr, stats::setNames(stats::plogis(stats::qlogis(params$survival) + bg), paste0("survival_", lev[2])))
  }
  if (!is.null(params$detection_by_occasion)) tr <- c(tr, stats::setNames(p, paste0("detection_", seq_len(T))))
  else {
    tr <- c(tr, detection = params$detection)
    if (has_cov) tr <- c(tr, detection_beta_effort = params$detection_beta_effort)
  }
  list(data = data, truth = tr, latent = list(alive = alive, entry = e, group = g, avail = avail, p = p))
}
