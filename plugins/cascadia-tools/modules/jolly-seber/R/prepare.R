# jolly-seber: turn data (from simulate() or read_data()) into model inputs.

js_bucket <- function(x, b = 50) as.integer(max(b, ceiling(x / b) * b))

js_prepare <- function(data, settings) {
  ch <- as.matrix(data$ch)
  storage.mode(ch) <- "integer"
  T <- ncol(ch); n <- nrow(ch)
  occ <- data$occasions
  if (is.null(occ)) occ <- data.frame(occasion = seq_len(T), label = as.character(seq_len(T)), interval = c(rep(1, T - 1), NA))
  dt <- as.numeric(occ$interval[seq_len(T - 1)])
  dt[!is.finite(dt)] <- 1
  notes <- character()
  # groups: code 0 = unknown, 1 = reference level, 2 = other level
  G <- 1L; gcode <- rep(0L, n); lev <- NULL
  if (!is.null(data$group) && isTRUE(settings$group_effect) && length(data$group_levels) == 2) {
    lev <- data$group_levels
    gcode <- match(data$group, lev); gcode[is.na(gcode)] <- 0L
    G <- 2L
    if (any(gcode == 0)) notes <- c(notes, sprintf("%d of %d captured animals have unknown %s; they are included as a mixture over %s.",
                                                    sum(gcode == 0), n, data$group_name, paste(lev, collapse = "/")))
  }
  # detection design matrix (levels first, then standardised slopes)
  det <- settings$detection
  if (identical(det, "auto")) det <- c("constant", data$detection_covariates)
  det <- unlist(det)
  scaling <- list()
  if ("time" %in% det) {
    Xp <- diag(T); colnames(Xp) <- paste0("detection_", seq_len(T)); Kp_lev <- T
  } else {
    covs <- setdiff(det, "constant")
    Xp <- matrix(1, T, 1, dimnames = list(NULL, "detection"))
    for (v in covs) {
      if (!v %in% names(occ)) stop(sprintf("detection covariate '%s' not found among occasion columns (%s)", v, paste(names(occ), collapse = ", ")))
      x <- as.numeric(occ[[v]])
      if (anyNA(x)) stop(sprintf("detection covariate '%s' has missing values", v))
      s <- stats::sd(x); if (!is.finite(s) || s == 0) stop(sprintf("detection covariate '%s' does not vary among occasions", v))
      scaling[[v]] <- list(center = mean(x), scale = s)
      Xp <- cbind(Xp, (x - mean(x)) / s); colnames(Xp)[ncol(Xp)] <- paste0("detection_beta_", v)
    }
    Kp_lev <- 1L
  }
  if (identical(settings$survival, "time")) {
    Xphi <- diag(T - 1); colnames(Xphi) <- paste0("survival_", seq_len(T - 1)); Kphi_lev <- T - 1
  } else {
    Xphi <- matrix(1, T - 1, 1, dimnames = list(NULL, "survival")); Kphi_lev <- 1L
  }
  M <- settings$M
  if (is.null(M)) M <- max(400L, 5L * n)
  M <- as.integer(M)
  if (M <= n) stop(sprintf("data augmentation M (%d) must exceed the number of animals captured (%d)", M, n))
  # unique (history, group) patterns + the all-zero pattern for M - n augmented rows
  key <- paste(apply(ch, 1, paste, collapse = ""), gcode)
  uk <- unique(key)
  first <- match(uk, key)
  pat <- rbind(ch[first, , drop = FALSE], rep(0L, T))
  pg <- c(gcode[first], 0L)
  w <- c(as.numeric(table(factor(key, levels = uk))), M - n)
  U <- nrow(pat)
  Upad <- js_bucket(U, 50)
  chp <- matrix(0, Upad, T); chp[seq_len(U), ] <- pat
  wp <- numeric(Upad); wp[seq_len(U)] <- w
  gp <- numeric(Upad); gp[seq_len(U)] <- pg
  list(ch = ch, n = n, T = T, dt = dt, occ = occ, G = G, gcode = gcode, group_levels = lev,
       group_name = data$group_name, Xp = Xp, Kp_lev = as.integer(Kp_lev), Xphi = Xphi,
       Kphi_lev = as.integer(Kphi_lev), scaling = scaling, M = M, pat = pat, pat_g = pg, pat_w = w,
       U = U, Upad = Upad, chp = chp, wp = wp, gp = gp, notes = notes,
       time_unit = if (is.null(data$time_unit)) "time unit" else data$time_unit,
       labels = as.character(occ$label))
}

# Augmentation size changes only the weight of the all-zero pattern.
js_set_M <- function(prep, M) {
  prep$M <- as.integer(M)
  prep$pat_w[prep$U] <- M - prep$n
  prep$wp[prep$U] <- M - prep$n
  prep
}
