# occu-single: data structure, standardisation and design matrices.
#
# The data list (what simulate()$data and read_data() both return):
#   model        "static" (MacKenzie et al. 2002) or "dynamic" (MacKenzie et al. 2003)
#   y            static: sites x surveys matrix; dynamic: sites x surveys x seasons array.
#                Values 0/1, NA = survey not done.
#   site_covs    data.frame, one row per site (may have zero columns)
#   obs_covs     named list of arrays with the same dims as y (survey-level covariates)
#   yearly_covs  named list of sites x seasons matrices (dynamic only: site-by-season covariates)
#   formulas     list(psi=, p=, col=, ext=) one-sided formula strings
#   site_ids     character, one per site; season_labels (dynamic)
#
# occ_prepare() turns that into numeric design arrays shared by every engine,
# so unmarked, Nimble and the naive control see identical covariates.

occ_default_formulas <- list(psi = "~ 1", p = "~ 1", col = "~ 1", ext = "~ 1")

occ_clean <- function(x) {
  x <- gsub("[^A-Za-z0-9_.]+", "_", x)
  gsub("^_+|_+$", "", x)
}

occ_dims <- function(data) {
  d <- dim(data$y)
  if (data$model == "static") c(M = d[1], J = d[2], T = 1L) else c(M = d[1], J = d[2], T = d[3])
}

# Standardise design-matrix columns (mean 0, SD 1 over the rows that enter
# the likelihood). Done after formula transforms, so log(x), I(x^2) and
# interactions work; 0/1 dummy columns and the intercept are left alone.
occ_colscale <- function(X, rows = rep(TRUE, nrow(X)), standardize = TRUE) {
  K <- ncol(X); sc <- matrix(c(0, 1), 2, K, dimnames = list(c("mean", "sd"), colnames(X)))
  if (standardize) for (k in seq_len(K)) {
    v <- X[rows, k]
    if (k == 1 || all(v %in% c(0, 1))) next
    s <- stats::sd(v); if (!is.finite(s) || s == 0) next
    sc[, k] <- c(mean(v), s)
  }
  sc
}
occ_apply_scale <- function(X, sc) sweep(sweep(X, 2, sc["mean", ], "-"), 2, sc["sd", ], "/")

occ_mm <- function(f, frame, what) {
  f <- stats::as.formula(f)
  if (length(f) != 2) stop("formula for ", what, " must be one-sided, e.g. '~ forest + elevation'")
  X <- stats::model.matrix(f, frame)
  if (!identical(colnames(X)[1], "(Intercept)"))
    stop("formula for ", what, " must include an intercept (do not use '0 +' or '- 1')")
  X
}

occ_prepare <- function(data, standardize = TRUE) {
  model <- match.arg(data$model, c("static", "dynamic"))
  dm <- occ_dims(data); M <- dm[["M"]]; J <- dm[["J"]]; T <- dm[["T"]]
  y <- array(as.numeric(data$y), c(M, J, T))
  w <- array(as.numeric(!is.na(y)), c(M, J, T))
  forms <- merge_list(occ_default_formulas, data$formulas)
  forms <- lapply(forms, function(f) if (inherits(f, "formula")) paste(deparse(f), collapse = "") else as.character(f))
  sc <- data$site_covs
  if (is.null(sc)) sc <- data.frame(row.names = seq_len(M))
  oc <- data$obs_covs; if (is.null(oc)) oc <- list()
  yc <- data$yearly_covs; if (is.null(yc)) yc <- list()
  oc <- lapply(oc, function(a) {
    if (is.factor(a) || is.character(a)) { lv <- if (is.factor(a)) levels(a) else sort(unique(stats::na.omit(as.vector(a))))
      structure(array(as.vector(as.character(a)), c(M, J, T)), levels = lv) }
    else array(as.numeric(a), c(M, J, T)) })
  yc <- lapply(yc, function(a) if (is.numeric(a)) matrix(as.numeric(a), M, T) else matrix(as.character(a), M, T))
  # observation covariates only matter where a survey was done
  for (k in names(oc)) oc[[k]][w == 0] <- NA

  scaling <- list(); raw <- list(); levels_ <- list()
  rawinfo <- function(k, v) { raw[[k]] <<- v; scaling[[k]] <<- c(mean = mean(v), sd = stats::sd(v)) }
  for (k in names(sc)) {
    if (is.character(sc[[k]])) sc[[k]] <- factor(sc[[k]])
    if (is.factor(sc[[k]])) { levels_[[k]] <- levels(sc[[k]]); next }
    rawinfo(k, sc[[k]])
  }
  for (k in names(yc)) {
    if (!is.numeric(yc[[k]])) { levels_[[k]] <- sort(unique(stats::na.omit(as.vector(yc[[k]])))); next }
    rawinfo(k, as.vector(yc[[k]][, seq_len(max(T - 1, 1))]))
  }
  for (k in names(oc)) {
    if (!is.numeric(oc[[k]])) { levels_[[k]] <- attr(oc[[k]], "levels"); next }
    rawinfo(k, oc[[k]][w == 1])
  }
  src <- list(site = names(sc), yearly = names(yc), obs = names(oc))
  allowed <- list(psi = src$site, col = c(src$site, src$yearly), ext = c(src$site, src$yearly),
                  p = c(src$site, src$yearly, src$obs))
  subs <- if (model == "static") c("psi", "p") else c("psi", "col", "ext", "p")
  for (s in subs) {
    v <- all.vars(stats::as.formula(forms[[s]]))
    bad <- setdiff(v, allowed[[s]])
    if (length(bad)) stop(sprintf("formula for %s uses %s, which is not available as a %s covariate (available: %s)",
                                  s, paste(bad, collapse = ", "),
                                  switch(s, psi = "site", p = "site/season/survey", "site/season"),
                                  paste(allowed[[s]], collapse = ", ")))
  }
  fill_na <- function(v, k) {  # values in unsurveyed cells are irrelevant; keep model.matrix from dropping rows
    if (is.numeric(v)) { v[is.na(v)] <- scaling[[k]][["mean"]]; v }
    else { v[is.na(v)] <- levels_[[k]][1]; factor(v, levels = levels_[[k]]) }
  }
  colscale <- list()
  # psi: one row per site
  Xpsi <- occ_mm(forms$psi, sc, "occupancy (psi)")
  colscale$psi <- occ_colscale(Xpsi, standardize = standardize); Xpsi <- occ_apply_scale(Xpsi, colscale$psi)
  # col/ext: rows (site, season t) for transitions t -> t+1, t = 1..T-1
  Xcol <- Xext <- NULL
  if (model == "dynamic") {
    g <- expand.grid(i = seq_len(M), t = seq_len(T - 1))
    fr <- sc[g$i, , drop = FALSE]
    for (k in names(yc)) fr[[k]] <- fill_na(yc[[k]][cbind(g$i, g$t)], k)
    Xc <- occ_mm(forms$col, fr, "colonisation"); Xe <- occ_mm(forms$ext, fr, "extinction")
    colscale$col <- occ_colscale(Xc, standardize = standardize); Xc <- occ_apply_scale(Xc, colscale$col)
    colscale$ext <- occ_colscale(Xe, standardize = standardize); Xe <- occ_apply_scale(Xe, colscale$ext)
    Xcol <- array(Xc, c(M, T - 1, ncol(Xc)), dimnames = list(NULL, NULL, colnames(Xc)))
    Xext <- array(Xe, c(M, T - 1, ncol(Xe)), dimnames = list(NULL, NULL, colnames(Xe)))
  }
  # p: rows (site, survey, season), site fastest
  g <- expand.grid(i = seq_len(M), j = seq_len(J), t = seq_len(T))
  fr <- sc[g$i, , drop = FALSE]
  for (k in names(yc)) fr[[k]] <- fill_na(yc[[k]][cbind(g$i, g$t)], k)
  for (k in names(oc)) fr[[k]] <- fill_na(oc[[k]][cbind(g$i, g$j, g$t)], k)
  Xp2 <- occ_mm(forms$p, fr, "detection (p)")
  surveyed <- as.vector(w) == 1
  colscale$p <- occ_colscale(Xp2, surveyed, standardize = standardize); Xp2 <- occ_apply_scale(Xp2, colscale$p)
  Xp2[!surveyed, -1] <- 0
  Xp <- array(Xp2, c(M, J, T, ncol(Xp2)), dimnames = list(NULL, NULL, NULL, colnames(Xp2)))
  cn <- list(psi = colnames(Xpsi), p = colnames(Xp2))
  if (model == "dynamic") { cn$col <- dimnames(Xcol)[[3]]; cn$ext <- dimnames(Xext)[[3]] }
  des <- list(model = model, M = M, J = J, T = T, y = y, w = w,
              Xpsi = Xpsi, Xp = Xp, Xcol = Xcol, Xext = Xext, colnames = cn, forms = forms,
              scaling = scaling, colscale = colscale, raw = raw, levels = levels_, sources = src,
              site_ids = if (is.null(data$site_ids)) sprintf("site%03d", seq_len(M)) else data$site_ids,
              season_labels = if (is.null(data$season_labels)) as.character(seq_len(T)) else data$season_labels)
  des$subs <- subs
  des
}

# Public parameter names for the coefficients of one submodel.
occ_coef_public <- function(sub, cols) {
  lab <- c(psi = "psi", p = "p", col = "col", ext = "ext")[[sub]]
  ifelse(cols == "(Intercept)", paste0(lab, "_logit_intercept"), paste0(lab, "_beta_", occ_clean(cols)))
}

# Intercept on the probability scale: the value when every standardised
# model term is at its mean (factors at their reference level).
occ_intercept_names <- function(model) {
  if (model == "static") c(psi = "occupancy_prob", p = "detection_prob")
  else c(psi = "initial_occupancy_prob", col = "colonisation_prob", ext = "extinction_prob", p = "detection_prob")
}

# ---- linear predictors on probability scale, given coefficient vectors ----
occ_probs <- function(des, b) {
  lp <- function(X, beta) {
    d <- dim(X); K <- d[length(d)]
    v <- matrix(X, ncol = K) %*% beta
    array(stats::plogis(v), d[-length(d)])
  }
  out <- list(psi = as.vector(stats::plogis(des$Xpsi %*% b$psi)), p = lp(des$Xp, b$p))
  if (des$model == "dynamic") { out$col <- lp(des$Xcol, b$col); out$ext <- lp(des$Xext, b$ext)
    dim(out$col) <- dim(out$ext) <- c(des$M, des$T - 1) }
  out
}

# Per-site, per-season emission probabilities of the observed detections
# given occupied (e1) and unoccupied (e0).
occ_emissions <- function(des, p) {
  y <- des$y; w <- des$w
  y0 <- ifelse(w == 1, y, 0)
  l1 <- w * (y0 * log(p) + (1 - y0) * log1p(-p))
  e1 <- exp(apply(l1, c(1, 3), sum))
  det <- apply(y0, c(1, 3), sum) > 0
  e0 <- ifelse(det, 0, 1)
  list(e1 = matrix(e1, des$M, des$T), e0 = matrix(e0, des$M, des$T))
}

# Forward-backward: log-likelihood and P(z[i,t] = 1 | data) for every site/season.
# Hand-written, so it doubles as an independent check of unmarked's likelihood.
occ_forward_backward <- function(des, pr) {
  em <- occ_emissions(des, pr$p); M <- des$M; T <- des$T
  a0 <- matrix(0, M, T); a1 <- matrix(0, M, T); ll <- numeric(M)
  u0 <- (1 - pr$psi) * em$e0[, 1]; u1 <- pr$psi * em$e1[, 1]
  s <- u0 + u1; ll <- ll + log(s); a0[, 1] <- u0 / s; a1[, 1] <- u1 / s
  if (T > 1) for (t in 2:T) {
    g <- pr$col[, t - 1]; e <- pr$ext[, t - 1]
    u0 <- (a0[, t - 1] * (1 - g) + a1[, t - 1] * e) * em$e0[, t]
    u1 <- (a0[, t - 1] * g + a1[, t - 1] * (1 - e)) * em$e1[, t]
    s <- u0 + u1; ll <- ll + log(s); a0[, t] <- u0 / s; a1[, t] <- u1 / s
  }
  b0 <- matrix(1, M, T); b1 <- matrix(1, M, T)
  if (T > 1) for (t in (T - 1):1) {
    g <- pr$col[, t]; e <- pr$ext[, t]
    v0 <- (1 - g) * em$e0[, t + 1] * b0[, t + 1] + g * em$e1[, t + 1] * b1[, t + 1]
    v1 <- e * em$e0[, t + 1] * b0[, t + 1] + (1 - e) * em$e1[, t + 1] * b1[, t + 1]
    s <- v0 + v1; b0[, t] <- v0 / s; b1[, t] <- v1 / s
  }
  post <- a1 * b1 / (a0 * b0 + a1 * b1)
  list(loglik = sum(ll), post = post)
}

# Marginal (unconditional) occupancy per site and season implied by the model.
occ_marginal_psi <- function(des, pr) {
  ps <- matrix(0, des$M, des$T); ps[, 1] <- pr$psi
  if (des$T > 1) for (t in 2:des$T)
    ps[, t] <- ps[, t - 1] * (1 - pr$ext[, t - 1]) + (1 - ps[, t - 1]) * pr$col[, t - 1]
  ps
}

# Simulate a new detection array from the model with the same missing pattern.
occ_simulate_y <- function(des, pr) {
  M <- des$M; T <- des$T
  z <- matrix(0, M, T); z[, 1] <- stats::rbinom(M, 1, pr$psi)
  if (T > 1) for (t in 2:T)
    z[, t] <- stats::rbinom(M, 1, ifelse(z[, t - 1] == 1, 1 - pr$ext[, t - 1], pr$col[, t - 1]))
  zz <- aperm(array(z, c(M, T, des$J)), c(1, 3, 2))
  y <- array(stats::rbinom(M * des$J * T, 1, zz * pr$p), c(M, des$J, T))
  y[des$w == 0] <- NA
  list(y = y, z = z)
}

# Split a coefficient vector/draw (named as in nimble: beta_psi[1], ...) into submodels.
occ_split_coef <- function(v, des) {
  out <- list()
  for (s in des$subs) out[[s]] <- unname(v[paste0("beta_", s, "[", seq_along(des$colnames[[s]]), "]")])
  out
}
occ_coef_internal <- function(des) unlist(lapply(des$subs, function(s)
  paste0("beta_", s, "[", seq_along(des$colnames[[s]]), "]")))
