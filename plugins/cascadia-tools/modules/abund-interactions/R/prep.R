# abund-interactions: turn the data list (from simulate() or read_data()) into
# design matrices shared by every engine.
#
# Data list (what simulate()$data and read_data() both return):
#   yA, yS     M x J matrices: number of sub-occasions (e.g. camera days) in
#              occasion j at site i on which the dominant (A) / subordinate (S)
#              species was detected. NA = occasion not surveyed.
#   K          M x J matrix: number of sub-occasions surveyed (camera-days
#              operating). 0 or NA = not surveyed. K = 1 everywhere means
#              plain detection/non-detection occasions.
#   site_covs  data.frame (M rows) of site covariates
#   obs_covs   named list of M x J matrices (occasion covariates)
#   formulas   list(lambda, r, psi, p): one-sided formulas as strings
#   modifier   name of the site covariate that modifies the interaction
#              (habitat), or NULL for a constant interaction
#   species    c(dominant = "...", subordinate = "...")
#   site_ids   character(M); coords optional data.frame(x, y)

abi_terms <- function(f) {
  f <- stats::as.formula(f)
  attr(stats::terms(f), "term.labels")
}

abi_clean <- function(x) gsub("[^A-Za-z0-9_]+", "_", gsub("[()]", "", x))

abi_prepare <- function(data, standardize = TRUE) {
  yA <- as.matrix(data$yA); yS <- as.matrix(data$yS)
  M <- nrow(yA); J <- ncol(yA)
  K <- if (is.null(data$K)) matrix(1, M, J) else as.matrix(data$K)
  K[is.na(K)] <- 0
  K[is.na(yA) | is.na(yS)] <- 0
  yA[K == 0] <- 0; yS[K == 0] <- 0
  sc <- if (is.null(data$site_covs) || !ncol(data$site_covs)) data.frame(.dummy = rep(0, M)) else as.data.frame(data$site_covs)
  scaling <- list()
  for (k in names(sc)) {
    if (is.character(sc[[k]]) || is.logical(sc[[k]])) sc[[k]] <- factor(sc[[k]])
    if (is.numeric(sc[[k]]) && standardize && k != ".dummy") {
      ce <- mean(sc[[k]], na.rm = TRUE); s <- stats::sd(sc[[k]], na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      scaling[[k]] <- list(center = ce, scale = s)
      sc[[k]] <- (sc[[k]] - ce) / s
    }
  }
  oc <- data$obs_covs
  long <- sc[rep(seq_len(M), J), , drop = FALSE]
  for (k in names(oc)) {
    v <- oc[[k]]
    if (is.numeric(v)) {
      w <- v[K > 0]
      ce <- mean(w, na.rm = TRUE); s <- stats::sd(w, na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      if (!standardize) { ce <- 0; s <- 1 }
      scaling[[k]] <- list(center = ce, scale = s)
      v <- (v - ce) / s; v[is.na(v)] <- 0
      long[[k]] <- as.vector(v)
    } else {
      v <- as.vector(v); v[is.na(v)] <- names(which.max(table(v)))
      long[[k]] <- factor(v)
    }
  }
  fm <- data$formulas
  for (s in c("lambda", "r", "psi", "p")) if (is.null(fm[[s]])) fm[[s]] <- "~ 1"
  mod <- data$modifier
  notes <- character()
  if (!is.null(mod) && nzchar(mod)) {
    if (!mod %in% names(sc)) stop("modifier '", mod, "' is not a site covariate")
    if (!mod %in% abi_terms(fm$psi)) {
      fm$psi <- paste(fm$psi, "+", mod)
      notes <- c(notes, sprintf("The habitat modifier '%s' was added as a main effect on subordinate occupancy (an interaction needs its main effect).", mod))
    }
  } else mod <- NULL
  mm <- function(f, df) {
    X <- stats::model.matrix(stats::as.formula(f), df)
    colnames(X) <- abi_clean(colnames(X))
    X
  }
  Xl <- mm(fm$lambda, sc); Xpsi <- mm(fm$psi, sc)
  Xr2 <- mm(fm$r, long); Xp2 <- mm(fm$p, long)
  Xr <- array(Xr2, c(M, J, ncol(Xr2))); Xp <- array(Xp2, c(M, J, ncol(Xp2)))
  Mod <- if (is.null(mod)) NULL else {
    X <- stats::model.matrix(stats::as.formula(paste("~", mod)), sc)[, -1, drop = FALSE]
    colnames(X) <- abi_clean(colnames(X)); X
  }
  for (X in list(Xl, Xpsi, Xr2, Xp2)) if (anyNA(X)) stop("missing values in covariates used by the model")
  list(M = M, J = J, K = K, yA = yA, yS = yS, Xl = Xl, Xr = Xr, Xpsi = Xpsi, Xp = Xp, Mod = Mod,
       colnames = list(lambda = colnames(Xl), r = colnames(Xr2), psi = colnames(Xpsi), p = colnames(Xp2),
                       mod = if (is.null(Mod)) character() else colnames(Mod)),
       formulas = fm, modifier = mod, scaling = scaling, site_covs_std = sc, notes = notes,
       species = if (is.null(data$species)) c(dominant = "Dominant species", subordinate = "Subordinate species") else data$species,
       site_ids = if (is.null(data$site_ids)) sprintf("S%03d", seq_len(M)) else data$site_ids,
       coords = data$coords, site_covs_raw = data$site_covs)
}

# Internal coefficient names (Nimble node names), in a fixed order.
abi_internal <- function(des) {
  P <- lapply(des$colnames, length)
  c(paste0("beta_lam[", seq_len(P$lambda), "]"), paste0("beta_r[", seq_len(P$r), "]"),
    paste0("beta_psi[", seq_len(P$psi), "]"), "gamma0",
    if (P$mod > 0) paste0("gamma_mod[", seq_len(P$mod), "]"),
    paste0("beta_p[", seq_len(P$p), "]"))
}

# Public (biological) names for the same coefficients; intercepts are
# reported back-transformed (abundance, probabilities).
abi_public <- function(des) {
  cn <- des$colnames
  c("dominant_mean_abundance", if (length(cn$lambda) > 1) paste0("lambda_beta_", cn$lambda[-1]),
    "dominant_detection_per_individual", if (length(cn$r) > 1) paste0("r_beta_", cn$r[-1]),
    "subordinate_occupancy_no_dominant", if (length(cn$psi) > 1) paste0("psi_beta_", cn$psi[-1]),
    "psi_beta_dominant_abundance",
    if (length(cn$mod)) paste0("psi_beta_dominant_abundance_x_", cn$mod),
    "subordinate_detection", if (length(cn$p) > 1) paste0("p_beta_", cn$p[-1]))
}

# Transform a draws matrix with internal names to public names.
abi_to_public <- function(des, m) {
  m <- as.matrix(m)
  int <- abi_internal(des); pub <- abi_public(des)
  out <- m[, int, drop = FALSE]
  colnames(out) <- pub
  out[, "dominant_mean_abundance"] <- exp(out[, "dominant_mean_abundance"])
  for (k in c("dominant_detection_per_individual", "subordinate_occupancy_no_dominant", "subordinate_detection"))
    out[, k] <- stats::plogis(out[, k])
  out
}

# Linear predictors from one coefficient vector (internal names)
abi_linpred <- function(des, th) {
  g <- function(pre, n) th[paste0(pre, "[", seq_len(n), "]")]
  P <- lapply(des$colnames, length)
  bl <- g("beta_lam", P$lambda); br <- g("beta_r", P$r); bpsi <- g("beta_psi", P$psi); bp <- g("beta_p", P$p)
  gm <- if (P$mod > 0) g("gamma_mod", P$mod) else numeric()
  lam <- exp(drop(des$Xl %*% bl))
  r <- stats::plogis(matrix(matrix(des$Xr, ncol = P$r) %*% br, des$M, des$J))
  p <- stats::plogis(matrix(matrix(des$Xp, ncol = P$p) %*% bp, des$M, des$J))
  eta <- drop(des$Xpsi %*% bpsi)
  slope <- th[["gamma0"]] + (if (P$mod > 0) drop(des$Mod %*% gm) else 0)
  list(lam = lam, r = r, p = p, eta = eta, slope = slope)
}
