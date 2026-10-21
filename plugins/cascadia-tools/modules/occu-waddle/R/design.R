# occu-waddle: parameter names, formula handling, design matrices, likelihood.
#
# Notation (Richmond et al. 2010; Waddle et al. 2010). Species A = dominant,
# species B = subordinate.
#   psiA   P(A occupies a site)
#   psiBA  P(B occupies | A present)          psiBa  P(B occupies | A absent)
#   pA     P(detect A on a survey | A present, B absent)
#   pB     P(detect B on a survey | B present, A absent)
#   rA     P(detect A | both present)
#   rBA    P(detect B | both present, A detected on that survey)
#   rBa    P(detect B | both present, A not detected on that survey)
# Each probability is logit-linear in covariates. A parameter without a
# formula is tied (wiqid / PRESENCE convention): psiBA -> psiBa, rA -> pA,
# rBa -> pB, rBA -> rBa.

wd_pars <- c("psiA", "psiBA", "psiBa", "pA", "pB", "rA", "rBA", "rBa")
wd_occ_pars <- c("psiA", "psiBA", "psiBa")
wd_tie <- c(psiBA = "psiBa", rA = "pA", rBa = "pB", rBA = "rBa")
wd_names <- c(psiA = "occ_dominant", psiBA = "occ_sub_given_dom_present", psiBa = "occ_sub_given_dom_absent",
              pA = "det_dominant_alone", pB = "det_sub_alone", rA = "det_dominant_with_sub",
              rBA = "det_sub_with_dom_detected", rBa = "det_sub_with_dom_missed")
wd_plain <- c(psiA = "occupancy of the dominant species",
              psiBA = "occupancy of the subordinate species where the dominant is present",
              psiBa = "occupancy of the subordinate species where the dominant is absent",
              pA = "detection of the dominant species where the subordinate is absent",
              pB = "detection of the subordinate species where the dominant is absent",
              rA = "detection of the dominant species where both occur",
              rBA = "detection of the subordinate species where both occur, on surveys that detected the dominant",
              rBa = "detection of the subordinate species where both occur, on surveys that missed the dominant")
wd_derived_names <- c("occ_subordinate", "occ_both", "species_interaction_factor",
                      "occ_log_odds_ratio_sub", "det_log_odds_ratio_sub")

wd_form_str <- function(f) {
  if (is.null(f)) return(NULL)
  f <- paste(trimws(as.character(f)), collapse = " ")
  if (!nzchar(f) || tolower(f) %in% c("tied", "same", "none", "na")) return(NULL)
  if (!startsWith(f, "~")) f <- paste("~", f)
  f
}

# Resolve formulas and ties. Returns list(forms = named formula strings of the
# free parameters, map = for each of the 8 parameters, the free parameter it
# uses).
wd_resolve <- function(formulas) {
  fl <- lapply(stats::setNames(wd_pars, wd_pars), function(p) wd_form_str(formulas[[p]]))
  for (p in c("psiA", "psiBa", "pA", "pB")) if (is.null(fl[[p]])) fl[[p]] <- "~ 1"
  bad <- setdiff(names(formulas), c(wd_pars, ""))
  if (length(bad)) stop("unknown formula key(s): ", paste(bad, collapse = ", "),
                        ". Use psiA, psiBA, psiBa, pA, pB, rA, rBA, rBa.")
  map <- stats::setNames(wd_pars, wd_pars)
  for (p in c("psiBA", "rA", "rBa", "rBA")) if (is.null(fl[[p]])) map[p] <- map[wd_tie[[p]]]
  free <- wd_pars[map[wd_pars] == wd_pars]
  for (p in free) {
    tt <- stats::terms(stats::as.formula(fl[[p]]))
    if (attr(tt, "intercept") != 1) stop("formula for ", p, " must keep the intercept")
    if (length(attr(tt, "offset"))) stop("offsets are not supported (", p, ")")
  }
  list(forms = fl[free], map = map, free = free)
}

wd_clean <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# Model matrix with numeric (non-0/1) columns standardised by the observed
# rows. Returns X (rows x K), and the pieces needed to rebuild it for newdata.
wd_mm <- function(form, df, use_rows) {
  f <- stats::as.formula(form)
  vars <- all.vars(f)
  miss <- setdiff(vars, names(df))
  if (length(miss)) stop("formula ", form, " uses unknown covariate(s): ", paste(miss, collapse = ", "))
  for (v in vars) if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  mf <- stats::model.frame(f, df[use_rows, , drop = FALSE], na.action = stats::na.fail)
  tt <- stats::delete.response(stats::terms(mf))
  xlev <- stats::.getXlevels(tt, mf)
  X <- wd_mm_raw(tt, xlev, df)
  ctr <- rep(0, ncol(X)); scl <- rep(1, ncol(X))
  for (k in seq_len(ncol(X))[-1]) {
    v <- X[use_rows, k]
    if (!all(v %in% c(0, 1))) { ctr[k] <- mean(v); scl[k] <- stats::sd(v); if (!is.finite(scl[k]) || scl[k] == 0) scl[k] <- 1 }
  }
  X <- sweep(sweep(X, 2, ctr, "-"), 2, scl, "/")
  X[!use_rows, ] <- 0
  X[is.na(X)] <- 0
  list(X = X, terms = tt, xlev = xlev, center = ctr, scale = scl, cols = colnames(X), vars = vars)
}

wd_mm_raw <- function(tt, xlev, df) {
  for (v in names(xlev)) df[[v]] <- factor(df[[v]], levels = xlev[[v]])
  mf <- stats::model.frame(tt, df, na.action = stats::na.pass, xlev = xlev)
  X <- stats::model.matrix(tt, mf)
  if (nrow(X) != nrow(df)) stop("internal: model matrix rows")
  X
}

# Build a model matrix for new covariate rows (raw units), standardised as in fit.
wd_mm_new <- function(mm, nd) {
  X <- wd_mm_raw(mm$terms, mm$xlev, nd)
  sweep(sweep(X, 2, mm$center, "-"), 2, mm$scale, "/")
}

# Everything the engines need from a data object.
wd_prepare <- function(data) {
  yA <- as.matrix(data$yA); yB <- as.matrix(data$yB)
  M <- nrow(yA); J <- ncol(yA)
  obs <- !is.na(yA) & !is.na(yB)
  rs <- wd_resolve(data$formulas)
  sc <- data$site_covs
  if (is.null(sc) || !ncol(sc)) sc <- data.frame(.dummy = rep(0, M))
  # survey-level data frame, row r = (i-1)*J + j
  sv <- sc[rep(seq_len(M), each = J), , drop = FALSE]
  for (k in names(data$obs_covs)) sv[[k]] <- as.vector(t(data$obs_covs[[k]]))
  sv$.survey <- rep(seq_len(J), M)
  obs_v <- as.vector(t(obs))
  mms <- list()
  for (p in rs$free) {
    if (p %in% wd_occ_pars) mms[[p]] <- wd_mm(rs$forms[[p]], sc, rep(TRUE, M))
    else mms[[p]] <- wd_mm(rs$forms[[p]], sv, obs_v)
  }
  # coefficient names (descriptive)
  cn <- list()
  for (p in rs$free) {
    cols <- mms[[p]]$cols
    cn[[p]] <- c(paste0(wd_names[[p]], "_intercept"),
                 if (length(cols) > 1) paste0(wd_names[[p]], "_beta_", wd_clean(cols[-1])))
  }
  K <- vapply(mms, function(m) ncol(m$X), 1L)
  idx <- split(seq_len(sum(K)), rep(names(K), K))[rs$free]
  ya <- as.vector(t(yA)); yb <- as.vector(t(yB))
  ya[!obs_v] <- 0; yb[!obs_v] <- 0
  list(M = M, J = J, yA = yA, yB = yB, obs = obs, obs_v = obs_v, ya = ya, yb = yb,
       forms = rs$forms, map = rs$map, free = rs$free, mms = mms, K = K, idx = idx,
       coef_names = unlist(cn, use.names = FALSE), site_covs = sc, surv_df = sv,
       species = data$species)
}

# Linear predictors -> probabilities for all 8 parameters (tied ones copy).
# beta: numeric vector of all coefficients; Xs: list of model matrices per free param.
wd_probs <- function(des, beta, Xs = lapply(des$mms, `[[`, "X")) {
  out <- list()
  for (p in des$free) out[[p]] <- stats::plogis(drop(Xs[[p]] %*% beta[des$idx[[p]]]))
  for (p in wd_pars) out[[p]] <- out[[des$map[[p]]]]
  out
}

# Hand-written log-likelihood in R (used by the maximum-likelihood engine and
# as the independent check of the Nimble likelihood). Sums site likelihoods
# over the four latent states (both, A only, B only, neither).
wd_loglik <- function(des, beta, per_site = FALSE) {
  q <- wd_probs(des, beta)
  M <- des$M; J <- des$J
  ya <- des$ya; yb <- des$yb; o <- des$obs_v
  lg <- function(p, y) ifelse(o, y * log(p) + (1 - y) * log1p(-p), 0)
  pbb <- ifelse(ya == 1, q$rBA, q$rBa)
  l1 <- rowSums(matrix(lg(q$rA, ya) + lg(pbb, yb), M, J, byrow = TRUE))
  l2 <- rowSums(matrix(lg(q$pA, ya), M, J, byrow = TRUE))
  l3 <- rowSums(matrix(lg(q$pB, yb), M, J, byrow = TRUE))
  anyA <- rowSums(matrix(ya, M, J, byrow = TRUE)) > 0
  anyB <- rowSums(matrix(yb, M, J, byrow = TRUE)) > 0
  lik <- q$psiA * q$psiBA * exp(l1) +
    ifelse(anyB, 0, q$psiA * (1 - q$psiBA) * exp(l2)) +
    ifelse(anyA, 0, (1 - q$psiA) * q$psiBa * exp(l3)) +
    ifelse(anyA | anyB, 0, (1 - q$psiA) * (1 - q$psiBa))
  ll <- log(lik)
  if (per_site) ll else sum(ll)
}

# Derived quantities from probabilities (vectors or scalars).
wd_derive <- function(psiA, psiBA, psiBa, rBA, rBa) {
  psiB <- psiA * psiBA + (1 - psiA) * psiBa
  cbind(occ_subordinate = psiB, occ_both = psiA * psiBA,
        species_interaction_factor = psiBA / psiB,
        occ_log_odds_ratio_sub = stats::qlogis(psiBA) - stats::qlogis(psiBa),
        det_log_odds_ratio_sub = stats::qlogis(rBA) - stats::qlogis(rBa))
}

# Which derived quantities are free (not fixed by a tie)?
wd_derived_free <- function(des) {
  keep <- wd_derived_names
  if (des$map[["psiBA"]] == des$map[["psiBa"]]) keep <- setdiff(keep, c("species_interaction_factor", "occ_log_odds_ratio_sub"))
  if (des$map[["rBA"]] == des$map[["rBa"]]) keep <- setdiff(keep, "det_log_odds_ratio_sub")
  keep
}

# Draws matrix of coefficients (columns = des$coef_names) -> reported
# columns: probabilities at average covariates (intercepts), link-scale
# slopes, derived quantities at average covariates.
wd_report_cols <- function(des, B, occ_only = FALSE) {
  B <- as.matrix(B)
  pars <- if (occ_only) intersect(des$free, wd_occ_pars) else des$free
  out <- list()
  for (p in pars) {
    ii <- des$idx[[p]]
    out[[wd_names[[p]]]] <- stats::plogis(B[, ii[1]])
    if (length(ii) > 1) for (k in ii[-1]) out[[des$coef_names[k]]] <- B[, k]
  }
  pr <- function(p) stats::plogis(B[, des$idx[[des$map[[p]]]][1]])
  if (occ_only) {
    d <- wd_derive(pr("psiA"), pr("psiBA"), pr("psiBa"), 0.5, 0.5)[, 1:4, drop = FALSE]
    keep <- intersect(colnames(d), wd_derived_free(des))
  } else {
    d <- wd_derive(pr("psiA"), pr("psiBA"), pr("psiBa"), pr("rBA"), pr("rBa"))
    keep <- wd_derived_free(des)
  }
  cbind(do.call(cbind, out), d[, keep, drop = FALSE])
}
