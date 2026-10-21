# occu-rota: data structure, latent-state bookkeeping and design matrices.
#
# The data list (what simulate()$data and read_data() both return):
#   species    character, S = 2 or 3 species names (letters, digits, _ only)
#   y          sites x surveys x species array of 0/1; NA = survey not done.
#              A survey is done for all species or none (one camera/visit
#              records every species), so the NA pattern is shared.
#   site_covs  data.frame, one row per site (may have zero columns)
#   obs_covs   named list of sites x surveys matrices (survey-level covariates)
#   formulas   list(occupancy = per-species one-sided formulas,
#                   interaction = per-pair (and triple) formulas, "0" = fixed at 0,
#                   detection = per-species formulas). Each element may be a
#                   single string (applied to every species / pair) or a named
#                   list keyed by species ("coyote") or pair ("coyote:fox").
#   max_order  highest interaction order (2 = pairwise; 3 allowed when S = 3)
#   site_ids   character, one per site
#
# rota_prepare() turns that into numeric design arrays shared by every
# engine, so unmarked, Nimble and the naive control see identical covariates.
#
# Natural-parameter model (Rota et al. 2016): latent state Z_i in {0,1}^S,
#   P(Z_i = z) = exp( sum_k f_ik * prod_{s in set_k} z_s ) / normaliser,
# with one natural parameter f_k per non-empty species set k up to max_order
# and f_ik = X_ik' beta_k. For two species: f1 = log(psi10/psi00),
# f2 = log(psi01/psi00), f12 = log(psi11 psi00 / (psi10 psi01)).

rota_clean <- function(x) {
  x <- gsub("[^A-Za-z0-9_.]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# All 2^S latent states: C x S 0/1 matrix, row 1 = no species present.
rota_states <- function(species) {
  S <- length(species)
  Z <- as.matrix(expand.grid(rep(list(0:1), S)))
  dimnames(Z) <- list(apply(Z, 1, function(r) if (!any(r == 1)) "none" else paste(species[r == 1], collapse = "+")),
                      species)
  storage.mode(Z) <- "double"
  Z
}

# Species sets (natural-parameter terms) up to max_order.
rota_terms <- function(species, max_order) {
  S <- length(species)
  out <- list()
  for (o in seq_len(min(max_order, S))) {
    cmb <- utils::combn(S, o, simplify = FALSE)
    for (cc in cmb) out[[paste(species[cc], collapse = ":")]] <- cc
  }
  out
}

# C x K matrix: A[c, k] = 1 if every species in term k is present in state c.
rota_Amat <- function(Z, terms) {
  A <- vapply(terms, function(cc) apply(Z[, cc, drop = FALSE], 1, prod), numeric(nrow(Z)))
  matrix(A, nrow(Z), length(terms), dimnames = list(rownames(Z), names(terms)))
}

# Natural parameters (n x K) -> state probabilities (n x C).
rota_psi <- function(Fm, A) {
  LW <- Fm %*% t(A)
  LW <- LW - apply(LW, 1, max)
  E <- exp(LW)
  P <- E / rowSums(E)
  colnames(P) <- rownames(A)
  P
}

# Biological-scale summaries of state probabilities (n x C), one column each:
# marginal occupancy, co-occurrence, conditional occupancy given another
# species present / absent (marginalising over any third species).
rota_derived_names <- function(species) {
  S <- length(species); nm <- paste0("occupancy_prob_", species)
  for (a in seq_len(S - 1)) for (b in (a + 1):S) nm <- c(nm, paste0("cooccurrence_prob_", species[a], "_", species[b]))
  for (a in seq_len(S)) for (b in seq_len(S)) if (a != b)
    nm <- c(nm, paste0("prob_", species[a], "_given_", species[b], c("_present", "_absent")))
  nm
}

rota_derived <- function(Psi, Z) {
  species <- colnames(Z); S <- length(species)
  Psi <- matrix(Psi, ncol = nrow(Z))
  out <- list()
  for (a in seq_len(S)) out[[paste0("occupancy_prob_", species[a])]] <- Psi %*% Z[, a]
  for (a in seq_len(S - 1)) for (b in (a + 1):S)
    out[[paste0("cooccurrence_prob_", species[a], "_", species[b])]] <- Psi %*% (Z[, a] * Z[, b])
  for (a in seq_len(S)) for (b in seq_len(S)) if (a != b) {
    pb <- Psi %*% Z[, b]
    out[[paste0("prob_", species[a], "_given_", species[b], "_present")]] <- (Psi %*% (Z[, a] * Z[, b])) / pb
    out[[paste0("prob_", species[a], "_given_", species[b], "_absent")]] <- (Psi %*% (Z[, a] * (1 - Z[, b]))) / (1 - pb)
  }
  m <- do.call(cbind, lapply(out, as.vector))
  m <- matrix(m, ncol = length(out)); colnames(m) <- names(out)
  m
}

# Standardise design-matrix columns (mean 0, SD 1 over the rows used), after
# formula transforms; intercept and 0/1 dummy columns are left alone.
rota_colscale <- function(X, rows = rep(TRUE, nrow(X)), standardize = TRUE) {
  K <- ncol(X); sc <- matrix(c(0, 1), 2, K, dimnames = list(c("mean", "sd"), colnames(X)))
  if (standardize) for (k in seq_len(K)) {
    if (colnames(X)[k] == "(Intercept)") next
    v <- X[rows, k]
    if (all(v %in% c(0, 1))) next
    s <- stats::sd(v); if (!is.finite(s) || s == 0) next
    sc[, k] <- c(mean(v), s)
  }
  sc
}
rota_apply_scale <- function(X, sc) sweep(sweep(X, 2, sc["mean", ], "-"), 2, sc["sd", ], "/")

rota_form_str <- function(f) if (inherits(f, "formula")) paste(deparse(f), collapse = "") else trimws(as.character(f))
rota_is_zero <- function(f) gsub("\\s", "", f) %in% c("0", "~0", "~-1", "~0+0")

rota_mm <- function(f, frame, what) {
  ff <- stats::as.formula(f)
  if (length(ff) != 2) stop("formula for ", what, " must be one-sided, e.g. '~ forest + elevation'")
  X <- stats::model.matrix(ff, frame)
  if (!identical(colnames(X)[1], "(Intercept)"))
    stop("formula for ", what, " must include an intercept (do not use '0 +' or '- 1'); use \"0\" to fix the whole term at zero")
  X
}

# Resolve the per-term formula: a single string applies to every species /
# pair; a named list is looked up by name (pairs in either order).
rota_pick_formula <- function(spec, key, default) {
  if (is.null(spec)) return(default)
  if (is.character(spec) && length(spec) == 1 && is.null(names(spec))) return(spec)
  spec <- as.list(spec)
  if (!is.null(spec[[key]])) return(rota_form_str(spec[[key]]))
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  for (nm in names(spec)) {
    p2 <- strsplit(nm, ":", fixed = TRUE)[[1]]
    if (length(p2) == length(parts) && setequal(trimws(p2), parts)) return(rota_form_str(spec[[nm]]))
  }
  if (!is.null(spec[["default"]])) return(rota_form_str(spec[["default"]]))
  default
}

rota_public_col <- function(cn) rota_clean(cn)

rota_prepare <- function(data, standardize = TRUE) {
  species <- as.character(data$species); S <- length(species)
  if (!S %in% 2:3) stop("occu-rota supports 2 or 3 species (got ", S, ")")
  y <- data$y
  M <- dim(y)[1]; J <- dim(y)[2]
  y <- array(as.numeric(y), c(M, J, S))
  w <- matrix(as.numeric(!is.na(y[, , 1])), M, J)
  max_order <- if (is.null(data$max_order)) 2L else as.integer(data$max_order)
  if (max_order < 2) stop("max_order must be 2 or 3")
  if (max_order > S) max_order <- S
  fm <- data$formulas; if (is.null(fm)) fm <- list()
  sc <- data$site_covs; if (is.null(sc)) sc <- data.frame(row.names = seq_len(M))
  for (k in names(sc)) if (is.character(sc[[k]])) sc[[k]] <- factor(sc[[k]])
  oc <- data$obs_covs; if (is.null(oc)) oc <- list()
  oc <- lapply(oc, function(a) { a <- matrix(a, M, J); a[w == 0] <- NA; a })

  levels_ <- list(); means <- list(); raw <- list()
  for (k in names(sc)) {
    if (is.factor(sc[[k]])) levels_[[k]] <- levels(sc[[k]])
    else { means[[k]] <- mean(sc[[k]]); raw[[k]] <- sc[[k]] }
  }
  for (k in names(oc)) {
    v <- as.vector(oc[[k]])
    if (is.numeric(v)) { means[[k]] <- mean(v, na.rm = TRUE); raw[[k]] <- v[!is.na(v)] }
    else levels_[[k]] <- sort(unique(stats::na.omit(v)))
  }

  # --- occupancy (natural parameters)
  terms <- rota_terms(species, max_order)
  forms <- list()
  for (k in names(terms)) {
    forms[[k]] <- if (length(terms[[k]]) == 1) rota_pick_formula(fm$occupancy, k, "~ 1")
                  else if (length(terms[[k]]) == 2) rota_pick_formula(fm$interaction, k, "~ 1")
                  else rota_pick_formula(fm$interaction, k, "0")
  }
  keep <- names(terms)[!vapply(forms, rota_is_zero, logical(1))]
  for (k in names(terms)[seq_len(S)]) if (!k %in% keep) stop("the occupancy formula for species ", k, " cannot be \"0\"")
  terms_all <- terms; terms <- terms[keep]; forms <- forms[keep]
  Xl <- list(); xscale <- list(); xcols <- list(); pub <- character(); term_of <- integer()
  for (ki in seq_along(terms)) {
    k <- names(terms)[ki]
    X <- rota_mm(forms[[k]], sc, paste0("'", k, "'"))
    s <- rota_colscale(X, standardize = standardize)
    X <- rota_apply_scale(X, s)
    Xl[[k]] <- X; xscale[[k]] <- s; xcols[[k]] <- colnames(X)
    sp <- strsplit(k, ":", fixed = TRUE)[[1]]
    stem <- if (length(sp) == 1) paste0("occ_", sp) else paste0("interaction_", paste(sp, collapse = "_"))
    pn <- ifelse(colnames(X) == "(Intercept)",
                 if (length(sp) == 1) paste0(stem, "_intercept") else stem,
                 paste0(stem, "_beta_", rota_public_col(colnames(X))))
    pub <- c(pub, pn); term_of <- c(term_of, rep(ki, ncol(X)))
  }
  Xall <- do.call(cbind, Xl); colnames(Xall) <- pub
  K <- length(terms); P <- ncol(Xall)

  # --- detection
  frame_det <- data.frame(row.names = seq_len(M * J))
  for (k in names(sc)) frame_det[[k]] <- rep(sc[[k]], each = J)        # site-major rows
  for (k in names(oc)) {
    v <- as.vector(t(oc[[k]]))
    frame_det[[k]] <- if (!is.null(levels_[[k]])) factor(v, levels = levels_[[k]]) else v
  }
  used <- as.vector(t(w)) == 1
  dforms <- list(); Wl <- list(); wscale <- list(); dpub <- character(); det_of <- integer()
  for (si in seq_len(S)) {
    s <- species[si]
    f <- rota_pick_formula(fm$detection, s, "~ 1")
    if (rota_is_zero(f)) stop("detection formula for ", s, " cannot be \"0\"")
    dforms[[s]] <- f
    fr <- frame_det
    vars <- all.vars(stats::as.formula(f))
    miss <- setdiff(vars, names(fr))
    if (length(miss)) stop("detection formula for ", s, " uses unknown covariate(s): ", paste(miss, collapse = ", "))
    # rows not surveyed: fill numeric NA with 0 / factor NA with first level so model.matrix keeps every row
    for (v in vars) {
      if (is.factor(fr[[v]])) fr[[v]][is.na(fr[[v]])] <- levels(fr[[v]])[1]
      else fr[[v]][is.na(fr[[v]])] <- 0
    }
    X <- rota_mm(f, fr, paste0("detection of ", s))
    sc_ <- rota_colscale(X, used, standardize = standardize)
    X <- rota_apply_scale(X, sc_)
    X[!used, ] <- 0
    Wl[[s]] <- X; wscale[[s]] <- sc_
    pn <- ifelse(colnames(X) == "(Intercept)", paste0("detection_prob_", s),
                 paste0("p_", s, "_beta_", rota_public_col(colnames(X))))
    dpub <- c(dpub, pn); det_of <- c(det_of, rep(si, ncol(X)))
  }
  Wall <- do.call(cbind, Wl); colnames(Wall) <- dpub
  Q <- ncol(Wall)
  W3 <- array(NA_real_, c(M, J, Q))
  for (q in seq_len(Q)) W3[, , q] <- matrix(Wall[, q], M, J, byrow = TRUE)

  Z <- rota_states(species)
  hold <- rota_hold_values(c(forms, dforms), means, raw)
  list(species = species, S = S, M = M, J = J, y = y, w = w, max_order = max_order,
       terms = terms, terms_all = terms_all, forms = forms, dforms = dforms, K = K, P = P, Q = Q,
       X = Xall, xlist = Xl, xscale = xscale, xcols = xcols, term_of = term_of,
       W = W3, wlist = Wl, wscale = wscale, det_of = det_of,
       Z = Z, A = rota_Amat(Z, terms), C = nrow(Z),
       levels = levels_, means = means, hold = hold, raw = raw, site_covs = sc, obs_covs = oc,
       site_ids = data$site_ids, standardize = standardize)
}

# "Average conditions": the raw value of each numeric covariate at which its
# model term equals the term's mean (so standardised terms are 0 and the
# intercepts describe that site). For a covariate used untransformed this is
# its arithmetic mean; used only as log(x) it is the geometric mean; any other
# single transform f(x) is solved for f(u) = mean(f(x)). If a covariate enters
# several different terms (e.g. x and I(x^2)) the arithmetic mean is used.
rota_hold_values <- function(forms, means, raw) {
  hold <- means
  for (v in names(means)) {
    exprs <- unique(unlist(lapply(forms, function(f) {
      if (rota_is_zero(f)) return(character())
      tl <- attr(stats::terms(stats::as.formula(f)), "term.labels")
      tl[vapply(tl, function(t) identical(all.vars(str2lang(t)), v), logical(1))]
    })))
    if (length(exprs) != 1 || identical(exprs, v)) next
    e <- str2lang(exprs); x <- raw[[v]]
    g <- function(u) eval(e, stats::setNames(list(u), v), baseenv())
    target <- tryCatch(mean(g(x)), error = function(err) NA_real_)
    if (!is.finite(target)) next
    sol <- tryCatch(stats::uniroot(function(u) g(u) - target, range(x), tol = 1e-10)$root, error = function(err) NULL)
    if (!is.null(sol)) hold[[v]] <- sol
  }
  hold
}

# Description of the average-conditions site, for reports.
rota_hold_text <- function(des, vars) {
  vars <- vars[vars %in% names(des$hold)]
  fv <- names(des$levels)[names(des$levels) %in% vars | names(des$levels) %in% unlist(lapply(c(des$forms, des$dforms), function(f) all.vars(stats::as.formula(f))))]
  num <- vapply(vars, function(v) {
    h <- des$hold[[v]]; m <- des$means[[v]]
    sprintf("%s = %s%s", v, signif(h, 3), if (abs(h - m) > 1e-8 * max(1, abs(m))) {
      tr <- unlist(lapply(c(des$forms, des$dforms), function(f) attr(stats::terms(stats::as.formula(f)), "term.labels")))
      if (any(grepl(paste0("^log\\(\\s*", v, "\\s*\\)$"), tr))) " (geometric mean, because the model uses its logarithm)" else " (value at which its transformed term is at its mean)"
    } else " (mean)")
  }, "")
  fac <- vapply(fv, function(v) sprintf("%s = %s (reference level)", v, des$levels[[v]][1]), "")
  paste(c(num, fac), collapse = ", ")
}

# Natural parameters (n_draws x n_rows x K) for new covariate rows.
# newdata: data.frame in RAW units; variables absent are held at their
# average-conditions value (des$hold, numeric) or reference level (factor).
rota_newframe <- function(des, vars, newdata) {
  n <- max(1, nrow(newdata))
  nd <- data.frame(row.names = seq_len(n))
  for (v in vars) {
    if (!is.null(des$levels[[v]])) {
      val <- if (v %in% names(newdata)) as.character(newdata[[v]]) else rep(des$levels[[v]][1], n)
      nd[[v]] <- factor(val, levels = des$levels[[v]])
    } else nd[[v]] <- if (v %in% names(newdata)) newdata[[v]] else rep(if (!is.null(des$hold[[v]])) des$hold[[v]] else des$means[[v]], n)
  }
  nd
}

rota_Xnew <- function(des, newdata) {
  blocks <- lapply(names(des$terms), function(k) {
    f <- stats::as.formula(des$forms[[k]])
    nd <- rota_newframe(des, all.vars(f), newdata)
    rota_apply_scale(stats::model.matrix(f, nd), des$xscale[[k]])
  })
  X <- do.call(cbind, blocks); colnames(X) <- colnames(des$X); X
}

rota_Wnew <- function(des, s, newdata) {
  f <- stats::as.formula(des$dforms[[s]])
  nd <- rota_newframe(des, all.vars(f), newdata)
  rota_apply_scale(stats::model.matrix(f, nd), des$wscale[[s]])
}

# beta draws (n x P) and new rows Xn (r x P) -> list over rows of n x K natural params
rota_F_from <- function(des, B, Xn) {
  lapply(seq_len(nrow(Xn)), function(r) {
    vapply(seq_len(des$K), function(k) {
      cols <- which(des$term_of == k)
      as.vector(B[, cols, drop = FALSE] %*% Xn[r, cols])
    }, numeric(nrow(B)))
  })
}

# Site-level natural parameters (M x K) for one beta vector
rota_F_sites <- function(des, beta) {
  vapply(seq_len(des$K), function(k) {
    cols <- which(des$term_of == k)
    as.vector(des$X[, cols, drop = FALSE] %*% beta[cols])
  }, numeric(des$M))
}

# Detection probabilities M x J x S for one alpha vector
rota_p_sites <- function(des, alpha) {
  p <- array(0, c(des$M, des$J, des$S))
  for (s in seq_len(des$S)) {
    cols <- which(des$det_of == s)
    eta <- matrix(0, des$M, des$J)
    for (q in cols) eta <- eta + des$W[, , q] * alpha[q]
    p[, , s] <- stats::plogis(eta) * des$w
  }
  p
}

# Marginal log-likelihood (sum over latent states), used by the naive and
# check routines and by the hand-written cross-check.
rota_loglik <- function(des, beta, alpha) {
  Psi <- rota_psi(rota_F_sites(des, beta), des$A)
  p <- rota_p_sites(des, alpha)
  y <- des$y; y[is.na(y)] <- 0
  ll <- matrix(0, des$M, des$C)
  for (c in seq_len(des$C)) for (s in seq_len(des$S)) {
    if (des$Z[c, s] == 1) ll[, c] <- ll[, c] + rowSums(y[, , s] * log(pmax(p[, , s], 1e-300)) +
                                                      (1 - y[, , s]) * log1p(-p[, , s]))
    else ll[, c] <- ll[, c] + ifelse(rowSums(y[, , s]) > 0, -Inf, 0)
  }
  lp <- ll + log(Psi)
  mx <- apply(lp, 1, max)
  sum(mx + log(rowSums(exp(lp - mx))))
}
