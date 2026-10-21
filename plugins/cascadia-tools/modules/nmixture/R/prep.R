# nmixture: settings, covariate standardisation and design matrices.
#
# Internal data structure (returned by simulate() and read_data()):
#   model        "static" (closed, Royle 2004) or "dynamic" (open, Dail & Madsen 2011)
#   y            static: sites x visits matrix of counts; dynamic: sites x visits x seasons array
#                (NA = visit not done)
#   site_covs    data.frame, one row per site
#   obs_covs     named list of arrays shaped like y (visit-level covariates)
#   yearly_covs  named list of sites x seasons matrices (dynamic only)
#   formulas     list(lambda, p[, gamma, omega]) as one-sided formula strings
#   site_ids, season_labels, read_issues

`%||%` <- function(a, b) if (is.null(a)) b else a

nm_default_formulas <- list(lambda = "~ 1", p = "~ 1", gamma = "~ 1", omega = "~ 1")

nm_settings <- list(
  engine = "unmarked",      # unmarked (ML, default) | nimble (Bayesian, static only) | rlik (hand-written ML, cross-check)
  mixture = "P",            # P (Poisson) | NB (negative binomial) | ZIP (zero-inflated Poisson)
  dynamics = "constant",    # dynamic model only: constant | notrend (unmarked::pcountOpen)
  standardize = TRUE,       # centre/scale numeric covariates (effects are per 1 SD)
  K = NULL,                 # upper bound for abundance per site; NULL = chosen automatically and checked
  K_tail_tol = 1e-4,        # auto-K: refit with larger K until P(N > K) < tol at the largest fitted mean
  k_sensitivity = TRUE,     # refit at larger K and report the change in abundance estimates
  k_mult = c(1.5, 2),
  compare_mixtures = TRUE,  # also fit P / NB / ZIP and report AIC (static, unmarked/rlik)
  gof_nboot = 200,          # parametric bootstrap (unmarked) or posterior predictive (nimble) replicates
  n_draws = 1000,           # parameter draws used for derived abundance (likelihood engines)
  seed = 1,
  # nimble
  niter = 6000, nburnin = 1000, nchains = 3, thin = 5,
  prior_sd_lambda_int = 3, prior_sd_p_int = 1.5, prior_sd_beta = 2, prior_sd_log_size = 2
)

nm_fast <- merge_list(nm_settings, list(k_sensitivity = FALSE, compare_mixtures = FALSE, gof_nboot = 0,
                                        n_draws = 400, niter = 3000, nburnin = 1000, thin = 2))

nm_dims <- function(data) {
  d <- dim(data$y)
  c(M = d[1], J = d[2], T = if (length(d) == 3) d[3] else 1L)
}

nm_scale <- function(v) {
  cen <- mean(v, na.rm = TRUE); sc <- stats::sd(as.vector(v), na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- 1
  list(center = cen, scale = sc)
}

nm_clean <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# Model matrix that keeps NA rows (they are masked later).
nm_mm <- function(form, df) {
  f <- stats::as.formula(form)
  if (!length(all.vars(f))) return(matrix(1, nrow(df), 1, dimnames = list(NULL, "(Intercept)")))
  mf <- stats::model.frame(f, df, na.action = stats::na.pass)
  X <- stats::model.matrix(f, mf)
  if (nrow(X) != nrow(df)) stop("could not build design matrix for ", form)
  X[is.na(X)] <- 0
  X
}

nm_prepare <- function(data, settings) {
  dm <- nm_dims(data); M <- dm[["M"]]; J <- dm[["J"]]; T <- dm[["T"]]
  model <- data$model %||% "static"
  forms <- merge_list(nm_default_formulas, data$formulas)
  y3 <- array(as.numeric(data$y), c(M, J, T))
  scaling <- list()
  std <- isTRUE(settings$standardize)
  sc <- if (is.null(data$site_covs) || !ncol(data$site_covs)) data.frame(row.names = seq_len(M)) else data$site_covs
  sc <- as.data.frame(sc)
  for (k in names(sc)) {
    if (is.character(sc[[k]]) || is.logical(sc[[k]])) sc[[k]] <- factor(sc[[k]])
    else if (is.numeric(sc[[k]]) && std) { s <- nm_scale(sc[[k]]); scaling[[k]] <- s; sc[[k]] <- (sc[[k]] - s$center) / s$scale }
  }
  oc <- list()
  for (k in names(data$obs_covs)) {
    a <- data$obs_covs[[k]]
    a <- array(a, c(M, J, T))
    if (is.numeric(a) && std) { s <- nm_scale(a); scaling[[k]] <- s; a <- (a - s$center) / s$scale }
    oc[[k]] <- a
  }
  yc <- list()
  for (k in names(data$yearly_covs)) {
    a <- matrix(data$yearly_covs[[k]], M, T)
    if (is.numeric(a) && std) { s <- nm_scale(a[, seq_len(max(T - 1, 1))]); scaling[[k]] <- s; a <- (a - s$center) / s$scale }
    yc[[k]] <- a
  }
  # long data frames in unmarked order: site slowest, then season, then visit
  long_o <- function(a) as.vector(aperm(array(a, c(M, J, T)), c(2, 3, 1)))
  odf <- data.frame(row.names = seq_len(M * J * T))
  for (k in names(sc)) odf[[k]] <- sc[[k]][rep(seq_len(M), each = J * T)]
  for (k in names(yc)) odf[[k]] <- as.vector(t(yc[[k]]))[rep(seq_len(M * T), each = J)]
  for (k in names(oc)) {
    v <- long_o(oc[[k]])
    odf[[k]] <- if (is.character(v)) factor(v) else v
  }
  ydf <- data.frame(row.names = seq_len(M * T))
  for (k in names(sc)) ydf[[k]] <- sc[[k]][rep(seq_len(M), each = T)]
  for (k in names(yc)) ydf[[k]] <- as.vector(t(yc[[k]]))
  # mask: visit done AND all detection covariates present
  pv <- all.vars(stats::as.formula(forms$p))
  w <- !is.na(y3)
  for (k in intersect(pv, names(oc))) w <- w & !is.na(oc[[k]])
  y3[!w] <- NA
  Xlam <- nm_mm(forms$lambda, sc)
  Xp_long <- nm_mm(forms$p, odf)
  np <- ncol(Xp_long)
  Xp <- aperm(array(Xp_long, c(J, T, M, np)), c(3, 1, 2, 4))
  des <- list(model = model, M = M, J = J, T = T, y = y3, w = w, formulas = forms,
              Xlam = Xlam, Xp = Xp, sc = sc, odf = odf, ydf = ydf, oc = oc, yc = yc,
              scaling = scaling, site_ids = data$site_ids %||% sprintf("site%03d", seq_len(M)),
              season_labels = data$season_labels %||% as.character(seq_len(T)),
              colnames = list(lam = colnames(Xlam), p = colnames(Xp_long)))
  if (model == "dynamic") {
    Xg <- nm_mm(forms$gamma, ydf); Xo <- nm_mm(forms$omega, ydf)
    des$Xgam <- aperm(array(Xg, c(T, M, ncol(Xg))), c(2, 1, 3))   # M x T x k (row t drives t -> t+1)
    des$Xom <- aperm(array(Xo, c(T, M, ncol(Xo))), c(2, 1, 3))
    des$colnames$gam <- colnames(Xg); des$colnames$om <- colnames(Xo)
  }
  des
}

# Names of the internal parameter vector (link scale), in a fixed order.
nm_theta_names <- function(des, mixture, dynamics = "constant") {
  n <- c(paste0("lam[", seq_along(des$colnames$lam), "]"), paste0("p[", seq_along(des$colnames$p), "]"))
  if (des$model == "dynamic") {
    if (dynamics != "notrend") n <- c(n, paste0("gam[", seq_along(des$colnames$gam), "]"))
    n <- c(n, paste0("om[", seq_along(des$colnames$om), "]"))
  }
  if (mixture == "NB") n <- c(n, "log_size")
  if (mixture == "ZIP") n <- c(n, "logit_zi")
  n
}

# Internal theta -> per-site quantities.
nm_unpack <- function(th, des, mixture, dynamics = "constant") {
  g <- function(pre) th[grep(paste0("^", pre, "\\["), names(th))]
  bl <- g("lam"); bp <- g("p")
  lam <- as.vector(exp(des$Xlam %*% bl))
  np <- length(bp)
  eta <- matrix(des$Xp, ncol = np) %*% bp
  p <- array(stats::plogis(eta), c(des$M, des$J, des$T))
  p[!des$w] <- 0
  out <- list(lam = lam, p = p,
              size = if (mixture == "NB") exp(th[["log_size"]]) else NA_real_,
              zi = if (mixture == "ZIP") stats::plogis(th[["logit_zi"]]) else 0)
  if (des$model == "dynamic" && des$T > 1) {
    Tm1 <- des$T - 1
    bo <- g("om")
    out$om <- matrix(stats::plogis(matrix(des$Xom[, seq_len(Tm1), , drop = FALSE], ncol = length(bo)) %*% bo), des$M, Tm1)
    if (dynamics == "notrend") {
      out$gam <- matrix(lam, des$M, Tm1) * (1 - out$om)   # unmarked "notrend": gamma = lambda (1 - omega)
    } else {
      bg <- g("gam")
      out$gam <- matrix(exp(matrix(des$Xgam[, seq_len(Tm1), , drop = FALSE], ncol = length(bg)) %*% bg), des$M, Tm1)
    }
  }
  out
}
