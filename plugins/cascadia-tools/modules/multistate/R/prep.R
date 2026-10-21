# multistate: turn data (simulated or read from files) into model inputs.
#
# data$y is an individuals x occasions integer matrix: 0 = not seen,
# 1..S = seen in state s, S+1 = seen but state not recorded (multievent).
#
# unknown_handling (settings):
#   "model"  (default) multievent model: state assignment is a parameter
#   "drop"   treat "seen, state unknown" as not seen (Arnason-Schwarz on the rest)
#   "assign" recode "seen, state unknown" to settings$naive_state (NEGATIVE CONTROL:
#            what happens when uncertain records are forced into a default state)

ms_prepare <- function(data, settings) {
  y <- as.matrix(data$y); storage.mode(y) <- "integer"
  states <- as.character(data$states); S <- length(states); T <- ncol(y)
  unk <- S + 1L
  handling <- if (is.null(settings$unknown_handling)) "model" else settings$unknown_handling
  handling <- match.arg(handling, c("model", "drop", "assign"))
  n_unknown <- sum(y == unk)
  if (handling == "drop") y[y == unk] <- 0L
  if (handling == "assign") {
    ns <- settings$naive_state
    k <- if (is.null(ns)) 1L else if (is.character(ns)) match(ns, states) else as.integer(ns)
    if (is.na(k)) stop("settings$naive_state must be one of the state names: ", paste(states, collapse = ", "))
    y[y == unk] <- k
  }
  model <- if (handling == "model" && (identical(data$model, "multievent") || any(y == unk))) "multievent" else "AS"
  if (!is.null(settings$model) && settings$model %in% c("AS", "multievent") && handling == "model") model <- settings$model
  if (model == "AS" && any(y == unk))
    stop("the data contain 'seen, state unknown' records but model = AS; use model: multievent (or unknown_handling: drop)")
  seen <- rowSums(y > 0) > 0
  keep <- which(seen)
  y <- y[keep, , drop = FALSE]
  N <- nrow(y)
  f <- apply(y, 1, function(r) which(r > 0)[1])
  e1 <- y[cbind(seq_len(N), f)]
  # Covariates -------------------------------------------------------------
  fm <- data$formulas; if (is.null(fm)) fm <- list()
  covs <- list(survival = as.character(unlist(fm$survival)), transition = as.character(unlist(fm$transition)),
               detection = as.character(unlist(fm$detection)))
  ic <- data$ind_covs; oc <- data$occ_covs
  scaling <- list(); level <- list()
  get_cov <- function(k) {
    if (!is.null(ic) && k %in% names(ic)) {
      x <- as.numeric(ic[[k]])[keep]; type <- "individual"
    } else if (!is.null(oc) && k %in% names(oc)) {
      x <- as.numeric(oc[[k]]); type <- "occasion"
    } else stop(sprintf("covariate '%s' is not a column of the individual or occasion covariate tables", k))
    # Values that the model actually uses: every individual; for occasion
    # covariates, occasions 2..T for detection and 1..T-1 for survival and
    # transitions. A missing value elsewhere (e.g. effort at the first
    # occasion, or temperature after the last one) is harmless and is set to
    # the mean; a missing value that is used is an error.
    used <- rep(TRUE, length(x))
    if (type == "occasion") {
      used[] <- FALSE
      if (k %in% covs$detection) used[seq_len(T)[-1]] <- TRUE
      if (k %in% c(covs$survival, covs$transition)) used[seq_len(T - 1)] <- TRUE
    }
    if (anyNA(x[used])) {
      where <- if (type == "occasion") paste("occasion(s)", paste(if (!is.null(data$occasion_labels)) data$occasion_labels[used & is.na(x)] else which(used & is.na(x)), collapse = ", ")) else
        sprintf("%d animal(s)", sum(is.na(x)))
      stop(sprintf("covariate '%s' is missing (or not a number) for %s; fill it in or drop the covariate", k, where))
    }
    ctr <- mean(x, na.rm = TRUE); scl <- stats::sd(x, na.rm = TRUE); if (!is.finite(scl) || scl == 0) scl <- 1
    scaling[[k]] <<- c(center = ctr, scale = scl); level[[k]] <<- type
    z <- (x - ctr) / scl; z[is.na(z)] <- 0
    list(x = z, type = type)
  }
  allk <- unique(unlist(covs))
  cv <- stats::setNames(lapply(allk, get_cov), allk)
  for (k in covs$detection) if (cv[[k]]$type != "occasion")
    stop(sprintf("detection covariate '%s' is an individual covariate; detection covariates must be per occasion (effort, weather)", k))
  T1 <- T - 1
  gam_by_ind <- any(vapply(c(covs$survival, covs$transition), function(k) cv[[k]]$type == "individual", TRUE))
  # X arrays: by individual N x T1 x K, or by occasion T1 x K (dummy zero column if K = 0)
  mk <- function(ks, by_ind) {
    K <- length(ks); Kd <- max(K, 1)
    if (by_ind) {
      X <- array(0, c(N, T1, Kd))
      for (j in seq_along(ks)) {
        c0 <- cv[[ks[j]]]
        X[, , j] <- if (c0$type == "individual") matrix(c0$x, N, T1) else matrix(c0$x[1:T1], N, T1, byrow = TRUE)
      }
    } else {
      X <- matrix(0, T1, Kd)
      for (j in seq_along(ks)) X[, j] <- cv[[ks[j]]]$x[1:T1]
    }
    X
  }
  Xphi <- mk(covs$survival, gam_by_ind); Xpsi <- mk(covs$transition, gam_by_ind)
  Kp <- length(covs$detection); Xp <- matrix(0, T, max(Kp, 1))
  for (j in seq_along(covs$detection)) Xp[, j] <- cv[[covs$detection[j]]]$x
  # Transition pairs (from, to), r != s
  pairs <- which(matrix(TRUE, S, S) & !diag(S), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
  colnames(pairs) <- c("from", "to")
  list(y = y, N = N, S = S, T = T, states = states, f = f, e1 = e1, model = model, handling = handling,
       n_unknown = n_unknown, keep = keep, ids = if (!is.null(data$ids)) data$ids[keep] else as.character(keep),
       covs = covs, scaling = scaling, level = level, gam_by_ind = gam_by_ind,
       Xphi = Xphi, Xpsi = Xpsi, Xp = Xp, Kphi = length(covs$survival), Kpsi = length(covs$transition), Kp = Kp,
       pairs = pairs, raw_covs = cv, occasion_labels = data$occasion_labels)
}

# Public names of parameters ----------------------------------------------
ms_pair_names <- function(des) sprintf("%s_to_%s", des$states[des$pairs[, 1]], des$states[des$pairs[, 2]])
