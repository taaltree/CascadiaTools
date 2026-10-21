# cjs: data structure, covariate scaling and design matrices shared by all
# engines (marked, Nimble marginal, Nimble latent) and by simulate().
#
# Data structure (what read_data() and simulate()$data return):
#   list(ch        = n x T integer matrix of 0/1 capture histories (rows = individuals),
#        ids       = character n, individual ids,
#        occasions = character T, occasion labels (e.g. years),
#        intervals = numeric T-1, time between occasions (survival is per unit time),
#        ind       = data.frame, n rows: group column, age_class column, individual covariates,
#        occ       = data.frame, T rows: occasion (environmental) covariates; row t is used for
#                    the survival interval t -> t+1 and for detection at occasion t,
#        group     = name of the grouping column in `ind` (or NULL),
#        age_col   = name of the age-at-marking column in `ind` (values adult/juvenile) or NULL,
#        formulas  = list(phi = "~ ...", p = "~ ..."),
#        species, units = named list of covariate units (for text), occasion_unit)
#
# Variables available to formulas:
#   time       occasion as a factor (one parameter per interval / occasion)
#   Time       linear trend (0, 1, 2, ...)
#   <group>    the grouping factor (e.g. sex)
#   tsm        time since marking: "first" interval after marking vs "later" (transients; Pradel et al. 1997)
#   age_class  juvenile for the first interval after marking as a juvenile, adult afterwards (phi only)
#   td         trap dependence: 1 if caught at the previous occasion (p only; Pradel 1993, Gimenez et al. 2003)
#   individual covariates (numeric: standardised internally; character: factor)
#   occasion covariates (numeric: standardised internally)
#   (1 | time) in the phi formula: random time effect on logit survival (Nimble only)

cjs_clean <- function(x) gsub("_+$", "", gsub("^_+", "", gsub("[^A-Za-z0-9]+", "_", x)))

cjs_rhs <- function(f) { f <- trimws(f); if (!startsWith(f, "~")) f <- paste("~", f); f }

# Split "(1 | time)" off a phi formula string.
cjs_parse_formula <- function(f) {
  f <- cjs_rhs(f)
  re <- grepl("\\(\\s*1\\s*\\|\\s*time\\s*\\)", f)
  if (grepl("\\|", gsub("\\(\\s*1\\s*\\|\\s*time\\s*\\)", "", f)))
    stop("only the random effect '(1 | time)' is supported (in the survival formula)")
  g <- gsub("\\+?\\s*\\(\\s*1\\s*\\|\\s*time\\s*\\)", "", f)
  g <- trimws(gsub("~\\s*\\+", "~", g))
  if (g == "~" || g == "") g <- "~ 1"
  list(formula = g, random_time = re)
}

# Centre/scale for every numeric covariate (individual: over individuals,
# occasion: over occasions with values). Used identically by simulate() and fit().
# 0/1 indicators are left as they are (centre 0, scale 1).
cjs_scale_one <- function(v, type) {
  if (all(v[!is.na(v)] %in% c(0, 1))) return(list(center = 0, scale = 1, type = type, binary = TRUE))
  s <- stats::sd(v, na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
  list(center = mean(v, na.rm = TRUE), scale = s, type = type, binary = FALSE)
}

cjs_scaling <- function(data) {
  sc <- list()
  for (k in names(data$ind)) {
    v <- data$ind[[k]]
    if (is.numeric(v) && !identical(k, data$group)) {
      sc[[k]] <- cjs_scale_one(v, "individual")
    }
  }
  for (k in names(data$occ)) {
    v <- data$occ[[k]]
    if (is.numeric(v)) {
      sc[[k]] <- cjs_scale_one(v, "occasion")
    }
  }
  sc
}

cjs_first <- function(ch) apply(ch, 1, function(r) { w <- which(r > 0); if (length(w)) w[1] else NA_integer_ })
cjs_last <- function(ch) apply(ch, 1, function(r) { w <- which(r > 0); if (length(w)) max(w) else NA_integer_ })

# "_"-prefixed factor levels make model.matrix column names readable
# ("sex_male", "time_2019", "tsm_first").
cjs_fac <- function(x, levels) factor(paste0("_", x), levels = paste0("_", levels))

cjs_group_levels <- function(data) {
  if (is.null(data$group)) return(NULL)
  g <- data$ind[[data$group]]
  if (is.factor(g)) levels(g) else sort(unique(as.character(g)))
}

# Long design data: one row per individual x interval (phi) or occasion (p).
# Row order: individual-major, so individual i occupies a contiguous block.
cjs_design_data <- function(data, scaling = cjs_scaling(data)) {
  ch <- data$ch; n <- nrow(ch); T <- ncol(ch); Tm1 <- T - 1
  first <- cjs_first(ch)
  occ_lab <- data$occasions
  ind <- data$ind
  base_ind <- data.frame(row.names = seq_len(n))
  for (k in names(ind)) {
    v <- ind[[k]]
    if (identical(k, data$group)) {
      base_ind[[k]] <- cjs_fac(as.character(v), cjs_group_levels(data))
    } else if (identical(k, data$age_col)) {
      next
    } else if (is.numeric(v)) {
      s <- scaling[[k]]; base_ind[[k]] <- (v - s$center) / s$scale
    } else {
      lv <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      base_ind[[k]] <- cjs_fac(as.character(v), lv)
    }
  }
  occ <- data$occ
  occ_s <- data.frame(row.names = seq_len(T))
  for (k in names(occ)) {
    v <- occ[[k]]
    if (is.numeric(v)) { s <- scaling[[k]]; occ_s[[k]] <- (v - s$center) / s$scale }
    else occ_s[[k]] <- cjs_fac(as.character(v), sort(unique(as.character(v))))
  }
  juv <- if (!is.null(data$age_col)) tolower(as.character(ind[[data$age_col]])) %in% c("juvenile", "juv", "j", "young", "hy", "subadult") else rep(FALSE, n)
  # phi rows (i, t), t = 1..T-1
  ii <- rep(seq_len(n), each = Tm1); tt <- rep(seq_len(Tm1), n)
  dphi <- base_ind[ii, , drop = FALSE]
  dphi$time <- cjs_fac(occ_lab[tt], occ_lab[seq_len(Tm1)])
  dphi$Time <- tt - 1
  dphi$tsm <- cjs_fac(ifelse(tt == first[ii], "first", "later"), c("later", "first"))
  dphi$age_class <- cjs_fac(ifelse(tt == first[ii] & juv[ii], "juvenile", "adult"), c("adult", "juvenile"))
  for (k in names(occ_s)) dphi[[k]] <- occ_s[[k]][tt]
  dphi$.i <- ii; dphi$.t <- tt
  dphi$.pre <- as.numeric(tt < first[ii])                    # before first capture: survival fixed to 1 (padding)
  dphi$.dt <- data$intervals[tt]
  # p rows (i, t), t = 1..T (row t = 1 is padding only)
  ii2 <- rep(seq_len(n), each = T); tt2 <- rep(seq_len(T), n)
  dp <- base_ind[ii2, , drop = FALSE]
  tl <- pmax(tt2, 2)
  dp$time <- cjs_fac(occ_lab[tl], occ_lab[2:T])
  dp$Time <- tl - 2
  prev <- ch[cbind(ii2, pmax(tt2 - 1, 1))]
  dp$td <- as.numeric(tt2 > 1 & prev > 0)
  dp$tsm <- cjs_fac(ifelse(tl - 1 == first[ii2], "first", "later"), c("later", "first"))
  for (k in names(occ_s)) dp[[k]] <- occ_s[[k]][tl]
  dp$.i <- ii2; dp$.t <- tt2
  dp$.upto <- as.numeric(tt2 <= first[ii2])                  # at or before first capture: p fixed to 1 (padding)
  rownames(dphi) <- NULL; rownames(dp) <- NULL
  list(phi = dphi, p = dp, first = first, n = n, T = T, scaling = scaling)
}

cjs_special_vars <- list(phi = c("time", "Time", "tsm", "age_class"), p = c("time", "Time", "tsm", "td"))

cjs_check_vars <- function(form, process, data) {
  v <- all.vars(stats::as.formula(form))
  avail <- c(cjs_special_vars[[process]], data$group, setdiff(names(data$ind), data$age_col), names(data$occ))
  bad <- setdiff(v, avail)
  if (length(bad))
    stop(sprintf("%s formula '%s' uses %s, which is not in the data. Available: %s", process, form,
                 paste(bad, collapse = ", "), paste(unique(avail), collapse = ", ")))
  if (process == "phi" && "td" %in% v) stop("td (trap dependence) acts on detection: put it in the p formula")
  if (process == "p" && "age_class" %in% v) stop("age_class is only supported in the survival (phi) formula")
  if ("age_class" %in% v && is.null(data$age_col)) stop("age_class in the formula needs an age-at-marking column (data: age_col)")
  invisible(v)
}

# Model matrix with readable names; returns X plus what is needed to build new rows.
cjs_mm <- function(form, df) {
  f <- stats::as.formula(form)
  vars <- all.vars(f)
  mf <- stats::model.frame(f, df, na.action = stats::na.pass)
  X <- stats::model.matrix(f, mf)
  if (anyNA(X)) stop("missing covariate values in rows used by formula ", form)
  fv <- vars[vapply(vars, function(v) is.factor(df[[v]]), TRUE)]
  xlev <- lapply(df[fv], levels)
  colnames(X) <- cjs_clean(colnames(X)); colnames(X)[colnames(X) == "Intercept"] <- "(Intercept)"
  # drop columns for factor levels that never occur (would be unidentifiable)
  keep <- colSums(abs(X)) > 0
  list(X = X[, keep, drop = FALSE], terms = stats::delete.response(stats::terms(mf)), xlev = xlev, keep = keep,
       formula = form, vars = vars)
}

# New design rows for prediction (template = rows of the design df).
cjs_newX <- function(mmo, newdata) {
  mf <- stats::model.frame(mmo$terms, newdata, xlev = mmo$xlev, na.action = stats::na.pass)
  X <- stats::model.matrix(mmo$terms, mf)
  X[, mmo$keep, drop = FALSE]
}

# Everything an engine needs.
cjs_prepare <- function(data, settings = list()) {
  fo <- merge_list(data$formulas, settings$formulas)
  if (is.null(fo$phi)) fo$phi <- "~ 1"
  if (is.null(fo$p)) fo$p <- "~ 1"
  pf <- cjs_parse_formula(fo$phi)
  if (grepl("\\|", fo$p)) stop("random effects are only supported for survival: '(1 | time)' in the phi formula")
  fo$phi <- pf$formula; fo$p <- cjs_rhs(fo$p)
  random_time <- pf$random_time || isTRUE(settings$random_time)
  cjs_check_vars(fo$phi, "phi", data); cjs_check_vars(fo$p, "p", data)
  fo$phi_user <- fo$phi
  if (!identical(settings$tsm_interact, FALSE)) fo$phi <- cjs_cross_tsm(fo$phi)
  dd <- cjs_design_data(data)
  confounded <- "time" %in% all.vars(stats::as.formula(fo$phi)) && "time" %in% all.vars(stats::as.formula(fo$p))
  # Fully time-dependent phi(t) p(t): the last survival and last recapture probability are
  # not separately identifiable (only their product). Fix p at the last occasion to 1, so the
  # last 'survival' parameter estimates the product phi[T-1] * p[T]; every other parameter is
  # unchanged (the constrained and unconstrained maxima coincide on the ridge). `.full = TRUE`
  # keeps the unconstrained parametrisation (used only to compute true values in simulate()).
  fixed_last <- confounded && !isTRUE(settings$.full)
  dd$p$.fixT <- 0
  if (fixed_last) {
    T <- dd$T
    lv <- levels(dd$p$time)
    dd$p$time[dd$p$.t == T] <- lv[1]
    dd$p$time <- factor(as.character(dd$p$time), levels = lv[-length(lv)])
    dd$p$.fixT <- as.numeric(dd$p$.t == T)
  }
  mphi <- cjs_mm(fo$phi, dd$phi); mp <- cjs_mm(fo$p, dd$p)
  if (random_time && "time" %in% mphi$vars) stop("phi formula cannot contain both 'time' and '(1 | time)'")
  names_phi <- paste0("phi:", colnames(mphi$X)); names_p <- paste0("p:", colnames(mp$X))
  list(data = data, dd = dd, formulas = fo, random_time = random_time,
       mphi = mphi, mp = mp, Xphi = mphi$X, Xp = mp$X,
       theta_names = c(names_phi, names_p, if (random_time) "phi:sd_time"),
       n = dd$n, T = dd$T, first = dd$first, scaling = dd$scaling,
       confounded_last = confounded, fixed_last = fixed_last)
}

# Transients (Pradel et al. 1997): survival in the first interval after marking is
# (1 - tau) x resident survival, a MULTIPLICATIVE change on the probability scale. An
# additive logit model (tsm + x) forces covariate effects to be equal in the first and later
# intervals and attenuates them (validated: -16 to -25% bias). So tsm is crossed with the
# other survival terms ("~ tsm * (sex + x)"; MARK's phi(a2 * ...) model): main effects are
# then the RESIDENT effects. age_class is first-interval-only already, so it is not crossed.
# A formula that already has tsm interactions is left as written. settings tsm_interact: false
# turns this off.
cjs_cross_tsm <- function(form) {
  f <- stats::as.formula(form)
  lab <- attr(stats::terms(f), "term.labels")
  if (!"tsm" %in% lab || any(grepl(":", lab) & grepl("\\btsm\\b", lab))) return(form)
  oth <- setdiff(lab, c("tsm", "age_class"))
  if (!length(oth)) return(form)
  paste0("~ tsm * (", paste(oth, collapse = " + "), ")", if ("age_class" %in% lab) " + age_class" else "")
}
