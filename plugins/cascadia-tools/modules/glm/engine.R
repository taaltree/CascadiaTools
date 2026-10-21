# GLM / GLMM engine for the ecotools `glm` module.
#
# Pipeline: glm_prepare() (formula parsing, scaling, response/offset set-up)
#        -> glm_fit_tmb() or glm_fit_glm() (one internal representation `fo`)
#        -> glm_estimates() (contract table on link + biological scales)
#        -> glm_dharma() / glm_r2() (model checks)
# Everything here is family-generic; family-specific behaviour lives in the
# small lookup functions at the top.

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.atomic(a) && is.na(a))) b else a

GLM_FAMILIES <- c("gaussian", "binomial", "poisson", "nbinom2", "beta", "Gamma", "zipoisson", "zinb")
COUNT_FAMILIES <- c("poisson", "nbinom2", "zipoisson", "zinb")

canon_family <- function(f) {
  f0 <- tolower(gsub("[^A-Za-z0-9]", "", as.character(f %||% "")))
  map <- c(gaussian = "gaussian", normal = "gaussian", lm = "gaussian", linear = "gaussian",
           binomial = "binomial", bernoulli = "binomial", logistic = "binomial",
           poisson = "poisson", nbinom2 = "nbinom2", negbin = "nbinom2", negativebinomial = "nbinom2",
           nb = "nbinom2", nbinom = "nbinom2", beta = "beta", betafamily = "beta", gamma = "Gamma",
           zipoisson = "zipoisson", zip = "zipoisson", zeroinflatedpoisson = "zipoisson",
           zinb = "zinb", zinbinom2 = "zinb", zinbinom = "zinb", zeroinflatednegativebinomial = "zinb",
           zeroinflatednegbin = "zinb")
  if (!f0 %in% names(map))
    stop("family '", f, "' is not supported. Use one of: ", paste(GLM_FAMILIES, collapse = ", "), call. = FALSE)
  unname(map[f0])
}

default_link <- function(fam) switch(fam, gaussian = "identity", binomial = "logit", beta = "logit", "log")
is_zi <- function(fam) fam %in% c("zipoisson", "zinb")

family_object <- function(fam, link) {
  switch(fam,
    gaussian = stats::gaussian(link = link),
    binomial = stats::binomial(link = link),
    poisson = , zipoisson = stats::poisson(link = link),
    nbinom2 = , zinb = glmmTMB::nbinom2(link = link),
    beta = glmmTMB::beta_family(link = link),
    Gamma = stats::Gamma(link = link))
}

# Name of the dispersion estimand (NA = family has none)
disp_name <- function(fam) switch(fam, gaussian = "residual_sd", nbinom2 = , zinb = "nb_size_theta",
                                  beta = "beta_precision_phi", Gamma = "gamma_cv", NA_character_)
default_dispersion <- function(fam) switch(fam, gaussian = 1, nbinom2 = , zinb = 2, beta = 10, Gamma = 0.5, NA_real_)

# Prefix of the exponentiated-coefficient estimand
ratio_prefix <- function(fam, link) {
  if (link == "logit") return("odds_ratio")
  if (link == "log") return(if (fam %in% COUNT_FAMILIES) "rate_ratio" else "ratio_of_means")
  NA_character_
}
ratio_words <- function(fam, link) switch(ratio_prefix(fam, link) %||% "none",
  odds_ratio = "odds ratio", rate_ratio = "rate ratio", ratio_of_means = "ratio of means", "difference")

link_inv <- function(link) stats::make.link(link)$linkinv

fmt <- function(x, k = 3) {
  vapply(x, function(v) {
    if (!is.finite(v)) return(as.character(v))
    if (abs(v) >= 1000) return(formatC(round(v), format = "d", big.mark = ","))
    if (abs(v) >= 100) return(as.character(round(v)))
    format(signif(v, k), scientific = FALSE, drop0trailing = TRUE)
  }, "")
}
clean_name <- function(x) { x <- gsub("[^A-Za-z0-9_]+", "_", x); x <- gsub("_+", "_", x); gsub("^_|_$", "", x) }
fmt_unit <- function(u) clean_name(gsub(".", "p", format(u, scientific = FALSE, drop0trailing = TRUE), fixed = TRUE))

# ---- formula -------------------------------------------------------------------
parse_glm_formula <- function(fstr) {
  if (is.null(fstr) || !nzchar(fstr)) stop("formula is missing (e.g. count ~ forest + (1|site))", call. = FALSE)
  f <- tryCatch(stats::as.formula(fstr), error = function(e)
    stop("could not read formula '", fstr, "': ", conditionMessage(e), call. = FALSE))
  if (length(f) != 3) stop("formula needs a response on the left: response ~ predictors", call. = FALSE)
  if (!is.name(f[[2]]))
    stop("the left side of the formula must be a single column name. For successes/trials data give the ",
         "successes column as the response and the trials column under `trials:` in the config.", call. = FALSE)
  if (grepl("offset\\(", fstr))
    stop("put the effort/offset column under `offset:` in the config instead of offset() in the formula", call. = FALSE)
  bars <- reformulas::findbars(f)
  lhs <- vapply(bars, function(b) deparse1(b[[2]]), "")
  if (any(lhs != "1"))
    stop("only random intercepts such as (1|site) are supported; found (", paste(lhs[lhs != "1"], collapse = ", "),
         " | ...). Random slopes are not validated in this module.", call. = FALSE)
  groups <- vapply(bars, function(b) deparse1(b[[3]]), "")
  if (any(grepl("[:/]", groups)))
    stop("nested/crossed grouping written as (1|a/b) or (1|a:b) is not supported; create one column per grouping factor ",
         "and use (1|a) + (1|b)", call. = FALSE)
  fixed <- reformulas::nobars(f)
  if (!inherits(fixed, "formula")) fixed <- stats::as.formula(paste(deparse1(f[[2]]), "~ 1"))
  rhs_vars <- if (length(fixed) == 3) all.vars(fixed[[3]]) else character()
  list(formula = f, fixed = fixed, response = as.character(f[[2]]), groups = groups,
       vars = setdiff(rhs_vars, groups))
}

# ---- preparation ---------------------------------------------------------------
glm_prepare <- function(data, settings = list()) {
  sp <- data$spec
  df <- data$df
  pf <- parse_glm_formula(sp$formula)
  fam <- canon_family(settings$force_family %||% sp$family)
  link <- if (!is.null(settings$force_family)) default_link(fam) else (sp$link %||% default_link(fam))
  resp <- pf$response
  vars <- pf$vars
  ref <- sp$reference %||% list()
  for (v in vars) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
    if (is.factor(df[[v]])) {
      df[[v]] <- droplevels(df[[v]])
      if (!is.null(ref[[v]])) df[[v]] <- stats::relevel(df[[v]], ref = as.character(ref[[v]]))
    }
  }
  for (g in pf$groups) df[[g]] <- factor(df[[g]])
  num_vars <- vars[vapply(vars, function(v) is.numeric(df[[v]]), logical(1))]
  fac_vars <- setdiff(vars, num_vars)
  sc <- sp$scale %||% "all"
  sc <- if (identical(sc, "all") || isTRUE(sc)) num_vars else if (identical(sc, "none") || isFALSE(sc)) character() else intersect(unlist(sc), num_vars)
  raw <- df
  scaling <- list()
  for (v in sc) {
    m <- mean(df[[v]]); s <- stats::sd(df[[v]])
    if (!is.finite(s) || s == 0) stop("covariate '", v, "' has no variation", call. = FALSE)
    df[[v]] <- (df[[v]] - m) / s
    scaling[[v]] <- c(mean = m, sd = s)
  }
  trials <- if (fam == "binomial") sp$trials else NULL
  lhs <- if (!is.null(trials)) { df$.eco_fail <- df[[trials]] - df[[resp]]; sprintf("cbind(%s, .eco_fail)", resp) } else resp
  if (!is.null(sp[["offset"]])) df$.eco_offset <- log(df[[sp[["offset"]]]])
  rhs_full <- deparse1(pf$formula[[3]])
  rhs_fixed <- if (length(pf$fixed) == 3) deparse1(pf$fixed[[3]]) else "1"
  rhs <- if (isTRUE(settings$drop_random)) rhs_fixed else rhs_full
  f_fit <- stats::as.formula(paste(lhs, "~", rhs, if (!is.null(sp[["offset"]])) "+ offset(.eco_offset)"))
  environment(f_fit) <- globalenv()
  groups <- if (isTRUE(settings$drop_random)) character() else pf$groups
  # fixed-effect design (terms carry predvars so poly()/I() predict correctly)
  rhs_f <- stats::as.formula(paste("~", rhs_fixed)); environment(rhs_f) <- globalenv()
  mf <- stats::model.frame(rhs_f, df)
  tt <- stats::terms(mf)
  xlev <- lapply(df[fac_vars], levels)
  X <- stats::model.matrix(tt, mf)
  labs <- vapply(colnames(X), coef_label, "", xlev = xlev)
  # reference row: scaled vars at 0 (= mean), unscaled numeric at their mean, factors at reference level
  ref_row <- df[1, vars, drop = FALSE]
  for (v in num_vars) ref_row[[v]] <- if (v %in% sc) 0 else mean(df[[v]])
  for (v in fac_vars) ref_row[[v]] <- factor(xlev[[v]][1], levels = xlev[[v]])
  prep <- list(df = df, raw = raw, spec = sp, pf = pf, fam = fam, link = link,
               famobj = family_object(fam, link), zif = if (is_zi(fam)) stats::as.formula(sp$zi %||% "~1") else NULL,
               f_fit = f_fit, rhs_fixed = rhs_fixed, groups = groups, vars = vars, num_vars = num_vars,
               fac_vars = fac_vars, scaled = sc, scaling = scaling, xlev = xlev, tt = tt,
               X_cols = colnames(X), labels_coef = labs, ref_row = ref_row, response = resp, trials = trials,
               offset = sp[["offset"]], offset_per = as.numeric(sp$offset_per %||% 1))
  ref_raw <- ref_row
  for (v in sc) ref_raw[[v]] <- scaling[[v]][["mean"]]
  prep$ref_raw <- ref_raw
  prep$X_ref <- design_matrix(prep, ref_row)[1, ]
  prep
}

coef_label <- function(cn, xlev) {
  if (cn == "(Intercept)") return("intercept")
  parts <- strsplit(cn, ":", fixed = TRUE)[[1]]
  parts <- vapply(parts, function(p) {
    for (v in names(xlev)) if (startsWith(p, v) && substring(p, nchar(v) + 1) %in% xlev[[v]])
      return(paste0(v, "_", substring(p, nchar(v) + 1)))
    p
  }, "")
  clean_name(paste(parts, collapse = "_x_"))
}

design_matrix <- function(prep, newdata) {
  for (v in prep$fac_vars) newdata[[v]] <- factor(as.character(newdata[[v]]), levels = prep$xlev[[v]])
  mf <- stats::model.frame(prep$tt, newdata, xlev = prep$xlev)
  stats::model.matrix(prep$tt, mf)
}

# model-scale newdata from raw-scale values
to_model_scale <- function(prep, nd) {
  for (v in names(prep$scaling)) if (v %in% names(nd))
    nd[[v]] <- (nd[[v]] - prep$scaling[[v]][["mean"]]) / prep$scaling[[v]][["sd"]]
  nd
}

# ---- engines -------------------------------------------------------------------
glm_fit_tmb <- function(prep, REML = FALSE) {
  reml <- isTRUE(REML) && length(prep$groups) > 0
  run <- function(reml, start = NULL) {
    w <- character()
    args <- list(prep$f_fit, data = prep$df, family = prep$famobj, ziformula = prep$zif %||% ~0, REML = reml)
    if (!is.null(start)) args$start <- start
    mm <- withCallingHandlers(do.call(glmmTMB::glmmTMB, args),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    list(m = mm, w = w, ok = mm$fit$convergence == 0 && isTRUE(mm$sdr$pdHess))
  }
  # ML first. For REML, start the variance parameters at the ML solution: from glmmTMB's
  # default start the REML optimiser can stop at a spurious zero-SD boundary optimum
  # (seen for beta GLMMs in validation; restarting from ML found a better REML optimum).
  fit_ml <- run(FALSE)
  m <- fit_ml$m; warns <- fit_ml$w; reml_note <- NULL
  if (reml) {
    pm <- fit_ml$m$fit$par
    st <- list(theta = unname(pm[names(pm) == "theta"]))
    if (any(names(pm) == "betadisp")) st$betadisp <- unname(pm[names(pm) == "betadisp"])
    if (any(names(pm) == "betazi")) st$betazi <- unname(pm[names(pm) == "betazi"])
    fit_reml <- tryCatch(run(TRUE, st), error = function(e) NULL)
    if (!is.null(fit_reml) && (fit_reml$ok || !fit_ml$ok)) { m <- fit_reml$m; warns <- fit_reml$w }
    else { reml <- FALSE; reml_note <- "The REML fit did not converge, so the maximum-likelihood fit is reported (random-effect SDs may be biased low by roughly 10-20% with ~30 groups)." }
  }
  par <- m$fit$par; nm <- names(par)
  fe <- glmmTMB::fixef(m)
  # Non-Gaussian REML: glmmTMB integrates the fixed effects out as random effects and reports
  # their joint MODE, which for non-Gaussian likelihoods is not the marginal MLE. In independent
  # verification this biased the NB GLMM intercept by about -0.047 (site SD 0.7) and did not shrink
  # with more sites (35 and 140 sites alike), so mean_response_at_reference was biased ~ -4.6%.
  # Fix (the GLMM analogue of LMM REML = GLS at the REML variances): keep the REML variance
  # parameters (RE SDs, dispersion) and re-estimate the fixed effects (and zero-inflation
  # coefficients) by ML with those parameters held at their REML values.
  m_cond <- NULL
  if (reml && prep$fam != "gaussian") {
    st <- list(theta = unname(par[nm == "theta"]))
    mp <- list(theta = factor(rep(NA, length(st$theta))))
    if (any(nm == "betadisp")) { st$betadisp <- unname(par[nm == "betadisp"]); mp$betadisp <- factor(rep(NA, length(st$betadisp))) }
    if (any(nm == "betazi")) st$betazi <- unname(par[nm == "betazi"])
    st$beta <- unname(fe$cond)
    cf <- tryCatch({
      w <- character()
      mm <- withCallingHandlers(glmmTMB::glmmTMB(prep$f_fit, data = prep$df, family = prep$famobj,
                                                 ziformula = prep$zif %||% ~0, start = st, map = mp),
        warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
      if (mm$fit$convergence == 0 && isTRUE(mm$sdr$pdHess)) list(m = mm, w = w) else NULL
    }, error = function(e) NULL)
    if (!is.null(cf)) m_cond <- cf$m
    else warns <- c(warns, "Fixed effects could not be re-estimated by ML at the REML variance parameters; REML joint-mode fixed effects are reported (the intercept may be slightly biased).")
  }
  if (!is.null(m_cond)) {
    fe <- glmmTMB::fixef(m_cond)
    Vc <- suppressWarnings(as.matrix(stats::vcov(m_cond, full = TRUE)))
    kc <- seq_len(length(fe$cond) + length(fe$zi))   # cond then zi come first in glmmTMB order
  }
  beta <- fe$cond
  if (!identical(names(beta), prep$X_cols))
    stop("internal: glmmTMB coefficient order differs from the design matrix (", paste(names(beta), collapse = ","), ")")
  # full covariance of (cond, zi, disp, theta) in glmmTMB order; works for ML and REML
  V <- tryCatch(suppressWarnings(as.matrix(stats::vcov(m, full = TRUE))), error = function(e) NULL)
  nb <- length(beta)
  if (is.null(V)) V <- matrix(NA_real_, nb + length(fe$zi) + sum(nm %in% c("betadisp", "theta")),
                              nb + length(fe$zi) + sum(nm %in% c("betadisp", "theta")))
  cn <- colnames(V) %||% character(ncol(V))
  ib <- seq_len(nb)
  izi <- grep("^zi~", cn); id <- grep("^disp~", cn); ith <- grep("^theta_", cn)
  if (!is.null(m_cond)) {
    # fixed-effect / zi block from the conditional ML fit; variance-parameter block from REML;
    # cross-covariances set to 0 (as for GLS fixed effects at REML variances in an LMM)
    if (!identical(colnames(Vc)[kc], cn[c(ib, izi)])) stop("internal: conditional-fit covariance order differs")
    V[c(ib, izi), ] <- 0; V[, c(ib, izi)] <- 0
    V[c(ib, izi), c(ib, izi)] <- Vc[kc, kc]
  }
  Vb <- V[ib, ib, drop = FALSE]; dimnames(Vb) <- list(names(beta), names(beta))
  disp <- NULL
  pd <- par[nm == "betadisp"]
  if (length(pd) == 1 && length(id) == 1) {
    k <- if (prep$fam == "Gamma") -0.5 else 1   # Gamma: betadisp = log(shape), CV = shape^-1/2
    disp <- c(log_est = k * pd[[1]], log_se = abs(k) * sqrt(V[id, id]))
  }
  zi <- NULL
  if (length(fe$zi)) zi <- list(beta = fe$zi, V = V[izi, izi, drop = FALSE], idx = izi)
  re <- NULL
  pt <- par[nm == "theta"]
  if (length(pt)) {
    gn <- names(glmmTMB::VarCorr(m)$cond)
    re <- data.frame(group = gn, log_sd = unname(pt), log_se = sqrt(diag(V)[ith]))
  }
  conv <- m$fit$convergence == 0 && isTRUE(m$sdr$pdHess) && all(is.finite(diag(Vb)))
  list(engine = "glmmTMB", model = m, model_sim = m_cond %||% m, reml_ml_fixef = !is.null(m_cond), beta = beta, V = Vb, Vfull = V, ib = ib, disp = disp, zi = zi, re = re,
       converged = conv, warnings = unique(c(warns, reml_note)), message = m$fit$message, REML = reml,
       pdHess = isTRUE(m$sdr$pdHess), logLik = as.numeric(stats::logLik(m)), nobs = stats::nobs(m))
}

glm_fit_glm <- function(prep) {
  if (length(prep$groups)) stop("stats::glm cannot fit random effects")
  warns <- character()
  m <- withCallingHandlers(
    if (prep$fam == "nbinom2") MASS::glm.nb(prep$f_fit, data = prep$df, link = log)
    else stats::glm(prep$f_fit, family = prep$famobj, data = prep$df),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  beta <- stats::coef(m); Vb <- stats::vcov(m)
  rdf <- m$df.residual
  disp <- switch(prep$fam,
    gaussian = c(log_est = 0.5 * log(summary(m)$dispersion), log_se = 1 / sqrt(2 * rdf)),
    Gamma = c(log_est = 0.5 * log(summary(m)$dispersion), log_se = 1 / sqrt(2 * rdf)),
    nbinom2 = c(log_est = log(m$theta), log_se = m$SE.theta / m$theta), NULL)
  list(engine = "stats::glm", model = m, beta = beta, V = Vb, Vfull = Vb, ib = seq_along(beta), disp = disp,
       zi = NULL, re = NULL, converged = isTRUE(m$converged) && all(is.finite(diag(Vb))),
       warnings = unique(warns), message = "", pdHess = TRUE,
       logLik = as.numeric(stats::logLik(m)), nobs = stats::nobs(m))
}

# ---- estimates -----------------------------------------------------------------
per_text <- function(prep, cn) {
  lab <- function(v) prep$spec$labels[[v]] %||% v
  if (cn == "(Intercept)") return("at reference")
  if (cn %in% prep$scaled) return(sprintf("1 SD of %s (%s)", lab(cn), fmt(prep$scaling[[cn]][["sd"]])))
  if (cn %in% prep$num_vars) return(sprintf("1 unit of %s", lab(cn)))
  for (v in prep$fac_vars) if (startsWith(cn, v) && substring(cn, nchar(v) + 1) %in% prep$xlev[[v]])
    return(sprintf("%s = %s vs %s", lab(v), substring(cn, nchar(v) + 1), prep$xlev[[v]][1]))
  if (grepl(":", cn)) return("interaction (product of scaled terms)")
  "see term"
}

unit_value <- function(prep, v) {
  u <- prep$spec$units[[v]]
  if (is.list(u)) u <- u$value
  if (is.null(u)) NULL else as.numeric(u)
}
unit_text <- function(prep, v) {
  u <- unit_value(prep, v)
  lab <- prep$spec$labels[[v]] %||% v
  ut <- prep$spec$unit_text[[v]] %||% (if (is.list(prep$spec$units[[v]])) prep$spec$units[[v]]$text)
  if (!is.null(ut)) sprintf("%s %s of %s", fmt(u), ut, lab) else sprintf("%s units of %s", fmt(u), lab)
}

# Monte-Carlo draws of (beta, betazi) from the asymptotic normal of the MLE
coef_draws <- function(fo, n, seed) {
  idx <- c(fo$ib, if (!is.null(fo$zi)) fo$zi$idx)
  mu <- c(fo$beta, if (!is.null(fo$zi)) fo$zi$beta)
  S <- fo$Vfull[idx, idx, drop = FALSE]
  if (!all(is.finite(S))) return(NULL)
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  set.seed(seed)
  d <- MASS::mvrnorm(n, mu, S)
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  nb <- length(fo$beta)
  list(beta = d[, seq_len(nb), drop = FALSE],
       zi = if (!is.null(fo$zi)) d[, nb + seq_along(fo$zi$beta), drop = FALSE])
}

# Expected response (incl. zero inflation, at offset_per units of effort) for model-scale newdata
predict_mu <- function(prep, fo, nd, draws = NULL) {
  X <- design_matrix(prep, nd)
  off <- if (!is.null(prep$offset)) log(prep$offset_per) else 0
  inv <- link_inv(prep$link)
  zi_p <- function(b) {
    if (is.null(fo$zi)) return(0)
    Z <- if (length(fo$zi$beta) == 1) matrix(1, nrow(nd), 1) else
      stats::model.matrix(stats::delete.response(stats::terms(prep$zif)), nd)
    stats::plogis(Z %*% b)
  }
  if (is.null(draws)) return(as.vector(inv(X %*% fo$beta + off) * (1 - zi_p(fo$zi$beta))))
  eta <- X %*% t(draws$beta) + off
  mu <- inv(eta)
  if (!is.null(fo$zi)) mu <- mu * (1 - zi_p(t(draws$zi)))
  mu
}

glm_estimates <- function(fo, prep, settings) {
  z <- stats::qnorm(0.975)
  rows <- list()
  add <- function(param, est, se, lo, hi, scale, per = "") rows[[length(rows) + 1]] <<-
    data.frame(param = param, estimate = est, sd = se, lower = min(lo, hi), upper = max(lo, hi),
               scale = scale, per = per, stringsAsFactors = FALSE)
  rp <- ratio_prefix(prep$fam, prep$link)
  lab <- prep$labels_coef
  se_all <- sqrt(pmax(diag(fo$V), 0))
  for (j in seq_along(fo$beta)) {
    cn <- names(fo$beta)[j]; b <- fo$beta[[j]]; s <- se_all[[j]]; L <- lab[[cn]]
    add(paste0("beta_", L), b, s, b - z * s, b + z * s, sprintf("link (%s)", prep$link), per_text(prep, cn))
  }
  for (j in seq_along(fo$beta)) {
    cn <- names(fo$beta)[j]; if (cn == "(Intercept)") next
    b <- fo$beta[[j]]; s <- se_all[[j]]; L <- lab[[cn]]
    if (!is.na(rp)) add(paste0(rp, "_", L), exp(b), exp(b) * s, exp(b - z * s), exp(b + z * s),
                        ratio_words(prep$fam, prep$link), per_text(prep, cn))
    u <- if (cn %in% prep$num_vars) unit_value(prep, cn) else NULL
    if (!is.null(u)) {
      k <- if (cn %in% prep$scaled) u / prep$scaling[[cn]][["sd"]] else u
      bb <- b * k; ss <- s * k
      nm <- paste0(if (is.na(rp)) "beta" else rp, "_", L, "_per_", fmt_unit(u))
      if (is.na(rp)) add(nm, bb, ss, bb - z * ss, bb + z * ss, "difference in mean", unit_text(prep, cn))
      else add(nm, exp(bb), exp(bb) * ss, exp(bb - z * ss), exp(bb + z * ss), ratio_words(prep$fam, prep$link), unit_text(prep, cn))
    }
  }
  # expected response at the reference point
  x0 <- prep$X_ref
  off <- if (!is.null(prep$offset)) log(prep$offset_per) else 0
  inv <- link_inv(prep$link)
  ref_txt <- reference_text(prep)
  if (is.null(fo$zi)) {
    eta <- sum(x0 * fo$beta) + off; se <- sqrt(max(drop(t(x0) %*% fo$V %*% x0), 0))
    mu <- inv(eta)
    add("mean_response_at_reference", mu, abs(stats::make.link(prep$link)$mu.eta(eta)) * se,
        inv(eta - z * se), inv(eta + z * se), "response", ref_txt)
  } else {
    dr <- coef_draws(fo, settings$n_draws, settings$seed)
    if (!is.null(dr)) {
      mu <- predict_mu(prep, fo, prep$ref_row, dr)[1, ]
      add("mean_response_at_reference", predict_mu(prep, fo, prep$ref_row), stats::sd(mu),
          stats::quantile(mu, 0.025), stats::quantile(mu, 0.975), "response", ref_txt)
    } else add("mean_response_at_reference", predict_mu(prep, fo, prep$ref_row), NA, NA, NA, "response", ref_txt)
  }
  if (!is.null(fo$disp) && !is.na(disp_name(prep$fam))) {
    d <- fo$disp
    add(disp_name(prep$fam), exp(d[["log_est"]]), exp(d[["log_est"]]) * d[["log_se"]],
        exp(d[["log_est"]] - z * d[["log_se"]]), exp(d[["log_est"]] + z * d[["log_se"]]), "dispersion",
        switch(prep$fam, gaussian = "residual SD", Gamma = "coefficient of variation",
               beta = "precision (larger = less variable)", "NB size (larger = closer to Poisson)"))
  }
  if (!is.null(fo$zi)) {
    se <- sqrt(pmax(diag(fo$zi$V), 0)); bz <- fo$zi$beta
    if (length(bz) == 1) {
      add("zero_inflation_prob", stats::plogis(bz), stats::plogis(bz) * (1 - stats::plogis(bz)) * se,
          stats::plogis(bz - z * se), stats::plogis(bz + z * se), "probability", "share of excess (structural) zeros")
    } else for (j in seq_along(bz)) {
      add(paste0("zi_beta_", clean_name(sub("^\\(Intercept\\)$", "intercept", names(bz)[j]))), bz[[j]], se[[j]],
          bz[[j]] - z * se[[j]], bz[[j]] + z * se[[j]], "link (logit), zero-inflation", "")
    }
  }
  if (!is.null(fo$re)) for (i in seq_len(nrow(fo$re))) {
    r <- fo$re[i, ]
    add(paste0("re_sd_", clean_name(r$group)), exp(r$log_sd), exp(r$log_sd) * r$log_se,
        exp(r$log_sd - z * r$log_se), exp(r$log_sd + z * r$log_se), sprintf("link (%s) SD", prep$link),
        sprintf("between-%s SD", r$group))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

reference_text <- function(prep) {
  parts <- character()
  if (length(prep$num_vars)) parts <- c(parts, "continuous predictors at their mean")
  for (v in prep$fac_vars) parts <- c(parts, sprintf("%s = %s", prep$spec$labels[[v]] %||% v, prep$xlev[[v]][1]))
  if (!is.null(prep$offset)) parts <- c(parts, sprintf("per %s %s", fmt(prep$offset_per), prep$spec$labels[[prep$offset]] %||% prep$offset))
  if (length(prep$groups)) parts <- c(parts, "typical (average) group")
  if (!length(parts)) "" else paste(parts, collapse = "; ")
}

# ---- model checks --------------------------------------------------------------
glm_dharma <- function(fo, prep, settings) {
  out <- tryCatch({
    r <- DHARMa::simulateResiduals(fo$model_sim %||% fo$model, n = settings$n_sim_resid, seed = settings$seed, plot = FALSE)
    p <- c(uniformity = DHARMa::testUniformity(r, plot = FALSE)$p.value)
    dt <- DHARMa::testDispersion(r, plot = FALSE)
    p["dispersion"] <- dt$p.value
    p["outliers"] <- suppressWarnings(DHARMa::testOutliers(r, plot = FALSE)$p.value)
    if (prep$fam %in% COUNT_FAMILIES) p["zero_inflation"] <- DHARMa::testZeroInflation(r, plot = FALSE)$p.value
    list(resid = r$scaledResiduals, pred = r$fittedPredictedResponse, p = p,
         disp_ratio = unname(dt$statistic), nsim = settings$n_sim_resid)
  }, error = function(e) list(error = conditionMessage(e)))
  out
}

dharma_notes <- function(dh, prep) {
  if (is.null(dh)) return(character())
  if (!is.null(dh$error)) return(paste("DHARMa residual checks could not be run:", dh$error))
  p <- dh$p
  s <- sprintf("Residual checks (DHARMa, %d simulations): uniformity p = %s; dispersion ratio %s (p = %s); outliers p = %s%s.",
               dh$nsim, fmt(p[["uniformity"]], 2), fmt(dh$disp_ratio, 3), fmt(p[["dispersion"]], 2),
               fmt(p[["outliers"]], 2),
               if (!is.na(p["zero_inflation"])) sprintf("; zero-inflation p = %s", fmt(p[["zero_inflation"]], 2)) else "")
  warn <- character()
  if (p[["dispersion"]] < 0.05) {
    warn <- c(warn, if (dh$disp_ratio > 1) {
      if (prep$fam %in% c("poisson", "zipoisson")) "Residuals are OVERDISPERSED: refit with family nbinom2 (or zinb) before interpreting intervals, which are too narrow now."
      else if (prep$fam == "binomial" && is.null(prep$trials)) "Dispersion test flagged for 0/1 data; this usually signals a missing predictor or grouping factor rather than true overdispersion."
      else "Residuals are overdispersed relative to the model: intervals may be too narrow; consider a grouping factor (random intercept) or a more flexible family."
    } else "Residuals are underdispersed: intervals are likely conservative (too wide).")
  }
  if (p[["uniformity"]] < 0.05) warn <- c(warn, "The residual distribution departs from what the model expects (uniformity test): check the QQ plot, the family and link, and missing nonlinear terms.")
  if (!is.na(p["zero_inflation"]) && p[["zero_inflation"]] < 0.05)
    warn <- c(warn, "There are more (or fewer) zeros than the model predicts: consider a zero-inflated family (zipoisson/zinb) or nbinom2.")
  if (p[["outliers"]] < 0.05) warn <- c(warn, "More extreme residuals (outliers) than expected: check data entry for the most extreme points.")
  if (length(warn)) warn <- c(warn, "With large samples these tests can flag deviations too small to matter; judge their size from the residual plots.")
  c(s, warn)
}

glm_r2 <- function(fo, prep) {
  r <- tryCatch(suppressWarnings(suppressMessages(performance::r2(fo$model))), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  v <- unlist(lapply(r, function(x) as.numeric(x)[1]))
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  names(v) <- clean_name(names(v))
  v
}
