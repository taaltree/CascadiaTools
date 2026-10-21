# cjs: simulate capture histories from a CJS model with sex, an individual
# covariate (body mass), an occasion covariate (winter severity), optional
# transients, juveniles, trap dependence and random temporal variation.
#
# The true model is defined on the standardised covariates of the simulated
# sample (cjs_scaling), exactly as fit() standardises them, so the true
# coefficients are directly comparable to the estimates.

cjs_simulate <- function(params, design, seed) {
  set.seed(seed)
  P <- function(k, d = 0) if (is.null(params[[k]])) d else params[[k]]
  T <- design$n_occasions; Tm1 <- T - 1
  nnew <- design$n_marked_per_occasion
  first <- rep(seq_len(Tm1), each = nnew)
  n <- length(first)
  yrs <- as.character(design$first_year + seq_len(T) - 1)
  sex <- ifelse(stats::runif(n) < design$prop_male, "male", "female")
  mass <- round(design$mass_mean + design$mass_sd * stats::rnorm(n), 1)
  winter <- round(pmax(0, design$winter_mean + design$winter_sd * stats::rnorm(T)), 0)
  juv <- stats::runif(n) < (design$prop_juvenile %||% 0)
  ind <- data.frame(sex = factor(sex, levels = c("female", "male")), body_mass = mass)
  age_col <- NULL
  if ((design$prop_juvenile %||% 0) > 0) { ind$age_at_marking <- ifelse(juv, "juvenile", "adult"); age_col <- "age_at_marking" }
  occ <- data.frame(winter_severity = winter)
  data <- list(ch = NULL, ids = sprintf("B%04d", seq_len(n)), occasions = yrs, intervals = rep(1, Tm1),
               ind = ind, occ = occ, group = "sex", age_col = age_col,
               formulas = list(phi = design$phi_formula, p = design$p_formula),
               species = "simulated animals", occasion_unit = "year",
               units = list(body_mass = "g", winter_severity = "days with snow cover"))
  sc <- cjs_scaling(data)
  zm <- (mass - sc$body_mass$center) / sc$body_mass$scale
  zw <- (winter - sc$winter_severity$center) / sc$winter_severity$scale
  male <- as.numeric(sex == "male")
  eps <- stats::rnorm(Tm1, 0, P("temporal_sd"))
  epsp <- stats::rnorm(T, 0, P("p_temporal_sd")); epsp[1] <- 0
  transient <- stats::runif(n) < P("prop_transient")
  ch <- matrix(0L, n, T)
  for (i in seq_len(n)) {
    f <- first[i]; ch[i, f] <- 1L; alive <- TRUE
    for (t in f:Tm1) {
      lp <- stats::qlogis(P("survival", 0.6)) + P("phi_beta_sex_male") * male[i] + P("phi_beta_winter_severity") * zw[t] +
        P("phi_beta_body_mass") * zm[i] + P("phi_beta_age_class_juvenile") * (juv[i] && t == f) + eps[t]
      ph <- if (transient[i] && t == f) 0 else stats::plogis(lp)
      alive <- alive && stats::runif(1) < ph
      if (!alive) break
      lq <- stats::qlogis(P("detection_prob", 0.5)) + P("p_beta_sex_male") * male[i] + P("p_beta_td") * ch[i, t] + epsp[t + 1]
      ch[i, t + 1] <- as.integer(stats::runif(1) < stats::plogis(lq))
    }
  }
  data$ch <- ch
  list(data = data, truth = cjs_truth(data, params, eps, epsp))
}

# True values of every estimand the fit formula produces.
cjs_truth <- function(data, params, eps, epsp) {
  P <- function(k, d = 0) if (is.null(params[[k]])) d else params[[k]]
  des <- cjs_prepare(data, list(.full = TRUE))
  th <- stats::setNames(numeric(length(des$theta_names)), des$theta_names)
  occ <- data$occasions; T <- length(occ)
  for (nm in des$theta_names) {
    proc <- sub(":.*", "", nm); col <- sub("^[^:]+:", "", nm)
    if (col == "(Intercept)") {
      th[nm] <- if (proc == "phi") stats::qlogis(P("survival", 0.6)) else stats::qlogis(P("detection_prob", 0.5))
    } else if (col == "sd_time") {
      th[nm] <- P("temporal_sd")
    } else if (grepl("^time_", col) && !grepl("_", sub("^time_", "", col)) ) {
      lab <- sub("^time_", "", col)
      if (proc == "phi") th[nm] <- eps[match(lab, occ)] - eps[1] else th[nm] <- epsp[match(lab, occ)] - epsp[2]
    } else if (col == "tsm_first") {
      s <- P("survival", 0.6); th[nm] <- stats::qlogis(max((1 - P("prop_transient")) * s, 1e-6)) - stats::qlogis(s)
    } else {
      th[nm] <- P(paste0(proc, "_beta_", col))
    }
  }
  if ("time" %in% des$mphi$vars) th["phi:(Intercept)"] <- th["phi:(Intercept)"] + eps[1]
  if ("time" %in% des$mp$vars) th["p:(Intercept)"] <- th["p:(Intercept)"] + epsp[2]
  tg <- cjs_targets(des)
  tr <- cjs_eval(matrix(th, 1, dimnames = list(NULL, names(th))), tg)[1, ]
  tr[grepl("^proportion_transient", names(tr))] <- P("prop_transient")
  # tsm main effect and tsm interactions have no single true value (transients are a mixture)
  tr <- tr[!grepl("^phi_beta_tsm_first", names(tr))]
  tr
}
