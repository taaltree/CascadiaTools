# Engines for sdm-presence-only. All return
#   list(coef, vcov, loglik, converged, engine, draws (matrix, columns named
#        by parameter), notes)
# Parameter names:
#   log_density              log intensity (individuals / km2) at mean covariates  [needs survey data]
#   intensity_beta_<term>    habitat slope on log intensity per 1 SD (or SD^2 for _sq)
#   po_intercept             log record intensity at mean covariates (= log_density + bias intercept; not interpretable)
#   bias_beta_<cov>          slope of log record (sampling) intensity per 1 SD of the bias covariate
#   logit_detection          per-visit detection probability (logit), repeated-visit surveys

sdm_quad <- function(prep) {
  if (is.null(prep[["quad_idx"]])) list(q = seq_len(nrow(prep$X)), w = prep$A) else list(q = prep[["quad_idx"]], w = prep[["quad_w"]])
}

# ---- IPP by down-weighted Poisson regression (Renner et al. 2015) -----------
# Quadrature points = grid cell centres (weight = cell area), presence points
# get weight 1e-6 and response 1/weight. glm then maximises the IPP
# log-likelihood sum_pres log(lambda b) - integral(lambda b).
sdm_engine_ppm <- function(prep, settings) {
  qd <- sdm_quad(prep)
  wp <- 1e-6
  X <- rbind(prep$X[prep$pres_cell, , drop = FALSE], prep$X[qd$q, , drop = FALSE])
  Z <- rbind(prep$Z[prep$pres_cell, , drop = FALSE], prep$Z[qd$q, , drop = FALSE])
  w <- c(rep(wp, length(prep$pres_cell)), qd$w)
  yy <- c(rep(1 / wp, length(prep$pres_cell)), rep(0, length(qd$q)))
  M <- cbind(1, X, Z)
  colnames(M) <- c("po_intercept", paste0("intensity_beta_", colnames(prep$X)),
                   if (ncol(Z)) paste0("bias_beta_", colnames(prep$Z)))
  # explicit start + tight tolerance: from glm's default start (y + 0.1) the tiny
  # presence weights make the deviance criterion stop after one iteration
  st <- c(log(length(prep$pres_cell) / sum(qd$w)), rep(0, ncol(M) - 1))
  m <- suppressWarnings(stats::glm.fit(M, yy, weights = w, start = st, family = stats::poisson(),
                                       control = stats::glm.control(epsilon = 1e-12, maxit = 100)))
  V <- tryCatch(chol2inv(qr.R(m$qr)[seq_len(m$rank), seq_len(m$rank)]), error = function(e) NULL)
  cf <- m$coefficients; names(cf) <- colnames(M)
  if (!is.null(V)) { piv <- m$qr$pivot[seq_len(m$rank)]; V2 <- matrix(NA, ncol(M), ncol(M)); V2[piv, piv] <- V; V <- V2 }
  if (!is.null(V)) dimnames(V) <- list(colnames(M), colnames(M))
  ll <- sum(M[seq_along(prep$pres_cell), ] %*% cf) - sum(qd$w * exp(M[-seq_along(prep$pres_cell), ] %*% cf))
  list(coef = cf, vcov = V, loglik = ll, converged = isTRUE(m$converged) && m$rank == ncol(M),
       engine = "ppm (down-weighted Poisson regression)")
}

# ---- Integrated likelihood (Fithian et al. 2015; Koshkina et al. 2017) ------
sdm_par_names <- function(prep, use_po, use_survey) {
  c(if (use_survey) "log_density", paste0("intensity_beta_", colnames(prep$X)),
    if (use_po) c("po_intercept", if (ncol(prep$Z)) paste0("bias_beta_", colnames(prep$Z))),
    if (use_survey && isTRUE(prep$survey$estimate_p)) "logit_detection")
}

sdm_nll_factory <- function(prep, use_po, use_survey) {
  nm <- sdm_par_names(prep, use_po, use_survey)
  kx <- ncol(prep$X); kz <- ncol(prep$Z)
  qd <- sdm_quad(prep)
  Xq <- prep$X[qd$q, , drop = FALSE]; Zq <- prep$Z[qd$q, , drop = FALSE]
  Xp <- prep$X[prep$pres_cell, , drop = FALSE]; Zp <- prep$Z[prep$pres_cell, , drop = FALSE]
  sv <- prep$survey
  unpack <- function(th) {
    th <- stats::setNames(th, nm)
    list(a0 = if (use_survey) th[["log_density"]] else 0,
         b = th[paste0("intensity_beta_", colnames(prep$X))],
         c0 = if (use_po) th[["po_intercept"]] else 0,
         d = if (use_po && kz) th[paste0("bias_beta_", colnames(prep$Z))] else numeric(),
         lp = if (use_survey && isTRUE(sv$estimate_p)) th[["logit_detection"]] else Inf)
  }
  val_grad <- function(th) {
    u <- unpack(th); ll <- 0; g <- stats::setNames(numeric(length(nm)), nm)
    bn <- paste0("intensity_beta_", colnames(prep$X))
    if (use_po) {
      lin_q <- u$c0 + drop(Xq %*% u$b) + (if (kz) drop(Zq %*% u$d) else 0)
      lin_p <- u$c0 + drop(Xp %*% u$b) + (if (kz) drop(Zp %*% u$d) else 0)
      mu <- qd$w * exp(lin_q)
      ll <- ll + sum(lin_p) - sum(mu)
      g["po_intercept"] <- length(lin_p) - sum(mu)
      g[bn] <- g[bn] + colSums(Xp) - colSums(Xq * mu)
      if (kz) { dn <- paste0("bias_beta_", colnames(prep$Z)); g[dn] <- colSums(Zp) - colSums(Zq * mu) }
    }
    if (use_survey) {
      lam <- exp(u$a0 + drop(sv$X %*% u$b)); uu <- lam * sv$A
      psi <- -expm1(-uu); dpsi <- exp(-uu) * uu          # d psi / d linear predictor
      p <- stats::plogis(u$lp); y <- sv$ydet; J <- sv$J
      q <- if (is.finite(u$lp)) (1 - p)^J else 0
      det <- y > 0
      li <- numeric(length(y)); dl <- numeric(length(y)); dlp <- numeric(length(y))
      if (is.finite(u$lp)) {
        li[det] <- log(psi[det]) + y[det] * log(p) + (J[det] - y[det]) * log1p(-p)
        dlp[det] <- y[det] - J[det] * p
      } else li[det] <- log(psi[det])
      dl[det] <- dpsi[det] / psi[det]
      L0 <- 1 - psi[!det] * (1 - q[!det])
      li[!det] <- log(L0)
      dl[!det] <- -(1 - q[!det]) * dpsi[!det] / L0
      if (is.finite(u$lp)) dlp[!det] <- -psi[!det] * J[!det] * p * q[!det] / L0
      ll <- ll + sum(li)
      g["log_density"] <- sum(dl)
      g[bn] <- g[bn] + colSums(sv$X * dl)
      if (isTRUE(sv$estimate_p)) g["logit_detection"] <- sum(dlp)
    }
    list(ll = ll, grad = g)
  }
  list(names = nm, fn = function(th) { v <- val_grad(th)$ll; if (!is.finite(v)) 1e10 else -v },
       gr = function(th) -val_grad(th)$grad, unpack = unpack)
}

sdm_engine_ml <- function(prep, settings, use_po = TRUE, use_survey = !is.null(prep$survey)) {
  if (use_survey && is.null(prep$survey)) stop("engine needs survey data (data: survey)")
  f <- sdm_nll_factory(prep, use_po, use_survey)
  th0 <- stats::setNames(rep(0, length(f$names)), f$names)
  if (use_po) {
    s0 <- tryCatch(sdm_engine_ppm(prep, settings)$coef, error = function(e) NULL)
    if (!is.null(s0)) th0[intersect(names(s0), names(th0))] <- s0[intersect(names(s0), names(th0))]
    else th0["po_intercept"] <- log(prep$n_pres / sum(prep$A))
  }
  if (use_survey) {
    nv <- mean(prep$survey$ydet > 0); a <- mean(prep$survey$A)
    th0["log_density"] <- log(max(-log(1 - min(max(nv, 0.02), 0.95)), 1e-3) / a)
    if ("logit_detection" %in% names(th0)) th0["logit_detection"] <- 0.5
  }
  o <- stats::optim(th0, f$fn, f$gr, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  if (o$convergence != 0) o <- stats::optim(o$par, f$fn, f$gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  H <- stats::optimHess(o$par, f$fn, f$gr)
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (!is.null(V)) dimnames(V) <- list(f$names, f$names)
  ok <- o$convergence == 0 && !is.null(V) && all(diag(V) > 0)
  list(coef = o$par, vcov = V, loglik = -o$value, converged = ok,
       engine = if (use_po && use_survey) "integrated (presence-only + survey, joint likelihood)"
                else if (use_survey) "survey only (occupancy with cloglog link)" else "ppm (joint-likelihood code)")
}

# ---- maxnet (Phillips et al. 2017) -------------------------------------------
sdm_maxnet_formula <- function(prep, classes) {
  hab <- prep$habitat; t <- hab
  if (grepl("q", classes)) t <- c(t, paste0("I(", hab, "^2)"))
  if (grepl("h", classes)) t <- c(t, paste0("hinge(", hab, ")"))
  if (grepl("p", classes) && length(hab) > 1) { m <- outer(hab, hab, paste, sep = ":"); t <- c(t, m[lower.tri(m)]) }
  t <- c(t, prep$bias)   # bias covariates enter linearly only
  f <- stats::as.formula(paste("~", paste(t, collapse = " + "), "-1"))
  environment(f) <- asNamespace("maxnet")
  f
}

sdm_maxnet_df <- function(prep, rows) {
  as.data.frame(cbind(prep$X[rows, prep$habitat, drop = FALSE], prep$Z[rows, , drop = FALSE]))
}

sdm_maxnet_link <- function(betas, newdata) {
  obj <- list(betas = betas, alpha = 0, entropy = 0)
  drop(maxnet:::predict.maxnet(obj, newdata, clamp = FALSE, type = "link"))
}

sdm_engine_maxnet <- function(prep, settings) {
  eco_require(c("maxnet", "glmnet"))
  classes <- settings$maxnet_classes
  if (length(prep$quadratic) && !grepl("q", classes)) classes <- paste0(classes, "q")
  f <- sdm_maxnet_formula(prep, classes)
  qd <- sdm_quad(prep)
  bg <- sdm_maxnet_df(prep, qd$q)
  one <- function(pc) {
    d <- rbind(sdm_maxnet_df(prep, pc), bg)
    pp <- c(rep(1, length(pc)), rep(0, nrow(bg)))
    m <- maxnet::maxnet(pp, d, f = f, regmult = settings$maxnet_regmult, addsamplestobackground = FALSE)
    m$betas
  }
  b0 <- one(prep$pres_cell)
  B <- settings$maxnet_boot
  set.seed(settings$seed)
  boots <- if (B > 0) lapply(seq_len(B), function(i) tryCatch(one(sample(prep$pres_cell, replace = TRUE)), error = function(e) NULL)) else list()
  boots <- Filter(Negate(is.null), boots)
  feat <- unique(c(names(b0), unlist(lapply(boots, names))))
  bm <- t(vapply(boots, function(b) { v <- stats::setNames(rep(0, length(feat)), feat); v[names(b)] <- b; v }, numeric(length(feat))))
  if (!length(boots)) bm <- matrix(numeric(), 0, length(feat), dimnames = list(NULL, feat))
  colnames(bm) <- feat
  # map linear / quadratic features to the shared parameter names
  map <- c(stats::setNames(paste0("intensity_beta_", prep$habitat), prep$habitat),
           stats::setNames(paste0("intensity_beta_", prep$habitat, "_sq"), paste0("I(", prep$habitat, "^2)")),
           if (length(prep$bias)) stats::setNames(paste0("bias_beta_", prep$bias), prep$bias))
  cf <- stats::setNames(rep(0, length(map)), map); hit <- intersect(names(b0), names(map)); cf[map[hit]] <- b0[hit]
  keep <- map[names(map) %in% c(prep$habitat, if (grepl("q", classes)) paste0("I(", prep$habitat, "^2)"), prep$bias)]
  cf <- cf[keep]
  dr <- matrix(0, nrow(bm), length(keep), dimnames = list(NULL, keep))
  for (k in names(keep)) if (k %in% colnames(bm)) dr[, keep[[k]]] <- bm[, k]
  b0f <- stats::setNames(rep(0, length(feat)), feat); b0f[names(b0)] <- b0
  list(coef = cf, vcov = NULL, boot = dr, loglik = NA_real_, converged = length(boots) >= 0.9 * B,
       engine = sprintf("maxnet (features '%s', regmult %g)", classes, settings$maxnet_regmult),
       betas = b0f, boot_betas = bm, classes = classes)
}
