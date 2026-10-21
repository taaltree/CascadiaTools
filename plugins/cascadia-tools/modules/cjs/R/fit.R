# cjs: from engine output to the contract estimates table (biological-scale
# survival and detection, link-scale coefficients), goodness of fit (R2ucare).

cjs_settings <- list(engine = "marked", formulation = "marginal", niter = 6000, nburnin = 2000, nchains = 3,
                     thin = 1, prior_sd = 2, cache = TRUE, seed = 1, gof = TRUE, chat_adjust = FALSE,
                     n_mvn = 4000)
cjs_fast <- merge_list(cjs_settings, list(gof = FALSE, niter = 3000, nburnin = 1000, n_mvn = 2000))

cjs_lab <- function(f) sub("^_", "", as.character(f))

# Template prediction row: numeric covariates at their mean (0 after
# standardisation), Time at mid-study, factors at reference level.
cjs_template <- function(df, mmo) {
  r <- df[1, , drop = FALSE]
  for (v in mmo$vars) {
    if (is.factor(df[[v]])) r[[v]] <- factor(levels(df[[v]])[1], levels = levels(df[[v]]))
    else if (v == "Time") r[[v]] <- mean(range(df$Time))
    else if (v == "td") r[[v]] <- 0
    else r[[v]] <- 0
  }
  r
}

# Rows of factor-level combinations used for the survival/detection estimands.
cjs_grid <- function(des, process) {
  mmo <- if (process == "phi") des$mphi else des$mp
  df <- if (process == "phi") des$dd$phi else des$dd$p
  if (process == "p") df <- df[df$.t > 1, , drop = FALSE]
  tmpl <- cjs_template(df, mmo)
  fv <- mmo$vars[vapply(mmo$vars, function(v) is.factor(df[[v]]), TRUE)]
  lv <- lapply(fv, function(v) {
    l <- levels(df[[v]])
    if (v == "tsm") l <- "_later"
    l
  })
  names(lv) <- fv
  if ("td" %in% mmo$vars) lv$td <- c(0, 1)
  g <- if (length(lv)) expand.grid(lv, stringsAsFactors = FALSE) else data.frame(row.names = 1)
  rows <- tmpl[rep(1, nrow(g)), , drop = FALSE]
  nm <- character(nrow(g))
  # name order: group / individual factors, age class, time
  ord <- c(setdiff(names(lv), c("age_class", "time", "td", "tsm")), intersect(c("age_class", "time"), names(lv)), intersect("td", names(lv)))
  for (v in names(lv)) {
    if (v == "td") rows[[v]] <- g[[v]] else rows[[v]] <- factor(g[[v]], levels = levels(df[[v]]))
  }
  stem <- if (process == "phi") "survival" else "detection"
  for (k in seq_len(nrow(g))) {
    parts <- character()
    for (v in ord) {
      if (v == "td") { if (g[[v]][k] == 1) parts <- c(parts, "caught_previous") }
      else if (v != "tsm") parts <- c(parts, cjs_clean(cjs_lab(g[[v]][k])))
    }
    nm[k] <- paste(c(stem, parts), collapse = "_")
  }
  X <- cjs_newX(mmo, rows)
  list(X = X, names = nm, rows = rows, grid = g)
}

# Target estimands: list(name, type, ...). Types: lc_phi / lc_p (plogis of
# x'beta), coef, sd, ratio (1 - phi_first/phi_later), prod (phi_last * p_last).
cjs_targets <- function(des) {
  tg <- list()
  iphi <- grep("^phi:", des$theta_names); ip <- grep("^p:", des$theta_names)
  iphi <- setdiff(iphi, which(des$theta_names == "phi:sd_time"))
  gp <- cjs_grid(des, "phi"); gq <- cjs_grid(des, "p")
  occ <- des$data$occasions; T <- des$T
  last_phi <- paste0("_", occ[T - 1]); last_p <- paste0("_", occ[T])
  for (k in seq_along(gp$names)) {
    if (des$confounded_last && "time" %in% names(gp$grid) && gp$grid$time[k] == last_phi) {
      # p at the last occasion is fixed to 1: this survival parameter IS the product phi[T-1] * p[T]
      if (isTRUE(des$fixed_last) && !("age_class" %in% names(gp$grid) && gp$grid$age_class[k] != "_adult"))
        tg[[length(tg) + 1]] <- list(name = sub("^survival", "survival_x_detection", gp$names[k]), type = "lc",
                                     idx = iphi, x = gp$X[k, ], process = "phi")
      next
    }
    tg[[length(tg) + 1]] <- list(name = gp$names[k], type = "lc", idx = iphi, x = gp$X[k, ], process = "phi")
  }
  for (k in seq_along(gq$names)) {
    if (des$confounded_last && "time" %in% names(gq$grid) && gq$grid$time[k] == last_p) next
    tg[[length(tg) + 1]] <- list(name = gq$names[k], type = "lc", idx = ip, x = gq$X[k, ], process = "p")
  }
  if (des$confounded_last && !isTRUE(des$fixed_last)) {
    # Only the product phi_{T-1} * p_T is identifiable in a fully time-dependent model.
    kk <- which(gp$grid$time == last_phi)
    for (k in kk) {
      # matching detection row: same non-time factors
      oth <- setdiff(names(gp$grid), c("time", "tsm", "age_class"))
      cand <- which(gq$grid$time == last_p)
      for (v in intersect(oth, names(gq$grid))) cand <- cand[gq$grid[[v]][cand] == gp$grid[[v]][k]]
      if ("td" %in% names(gq$grid)) cand <- cand[gq$grid$td[cand] == 0]
      if ("age_class" %in% names(gp$grid) && gp$grid$age_class[k] != "_adult") next
      if (!length(cand)) next
      nm <- sub("^survival", "survival_x_detection", gp$names[k])
      tg[[length(tg) + 1]] <- list(name = nm, type = "prod", idx = iphi, x = gp$X[k, ], idx2 = ip, x2 = gq$X[cand[1], ])
    }
  }
  if ("tsm" %in% des$mphi$vars) {
    base <- which(!duplicated(gp$names))
    rows_l <- gp$rows; rows_f <- gp$rows
    rows_f$tsm <- factor("_first", levels = levels(des$dd$phi$tsm))
    Xf <- cjs_newX(des$mphi, rows_f)
    # one ratio per group level (other factors at the reference level)
    keep <- rep(TRUE, nrow(gp$grid))
    for (v in setdiff(names(gp$grid), c(des$data$group, "tsm"))) keep <- keep & gp$grid[[v]] == gp$grid[[v]][1]
    for (k in which(keep)) {
      suffix <- if (!is.null(des$data$group) && des$data$group %in% names(gp$grid)) paste0("_", cjs_clean(cjs_lab(gp$grid[[des$data$group]][k]))) else ""
      tg[[length(tg) + 1]] <- list(name = paste0("proportion_transient", suffix), type = "ratio", idx = iphi,
                                   x = gp$X[k, ], x2 = Xf[k, ])
    }
  }
  # Coefficients. Dummies of the `time` factor are not reported: survival_<interval> and
  # detection_<occasion> carry that information on the probability scale, and the logit
  # contrasts against the first interval are dominated by boundary years (and, with p fixed
  # at the last occasion, the last one is the product's coefficient).
  is_time <- function(nm) grepl("^(phi|p):(tsm_first_)?time_", nm)
  for (j in iphi) if (des$theta_names[j] != "phi:(Intercept)" && !is_time(des$theta_names[j]))
    tg[[length(tg) + 1]] <- list(name = paste0("phi_beta_", sub("^phi:", "", des$theta_names[j])), type = "coef", idx = j)
  for (j in ip) if (des$theta_names[j] != "p:(Intercept)" && !is_time(des$theta_names[j]))
    tg[[length(tg) + 1]] <- list(name = paste0("p_beta_", sub("^p:", "", des$theta_names[j])), type = "coef", idx = j)
  if (des$random_time)
    tg[[length(tg) + 1]] <- list(name = "survival_temporal_sd", type = "sd", idx = which(des$theta_names == "phi:sd_time"))
  nm <- vapply(tg, `[[`, "", "name")
  if (anyDuplicated(nm)) tg <- tg[!duplicated(nm)]
  tg
}

# Draws (rows of TH) -> matrix of target values.
cjs_eval <- function(TH, tg) {
  TH <- as.matrix(TH)
  out <- vapply(tg, function(g) {
    switch(g$type,
           lc = stats::plogis(drop(TH[, g$idx, drop = FALSE] %*% g$x)),
           coef = TH[, g$idx], sd = TH[, g$idx],
           ratio = 1 - stats::plogis(drop(TH[, g$idx, drop = FALSE] %*% g$x2)) / stats::plogis(drop(TH[, g$idx, drop = FALSE] %*% g$x)),
           prod = stats::plogis(drop(TH[, g$idx, drop = FALSE] %*% g$x)) * stats::plogis(drop(TH[, g$idx2, drop = FALSE] %*% g$x2)))
  }, numeric(nrow(TH)))
  out <- matrix(out, nrow(TH))
  colnames(out) <- vapply(tg, `[[`, "", "name")
  out
}

cjs_mvn <- function(theta, V, n) {
  ok <- all(is.finite(V))
  if (!ok) return(matrix(theta, n, length(theta), byrow = TRUE, dimnames = list(NULL, names(theta))))
  ev <- eigen(V, symmetric = TRUE)
  A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(theta))
  out <- sweep(matrix(stats::rnorm(n * length(theta)), n) %*% t(A), 2, theta, "+")
  colnames(out) <- names(theta); out
}

# Likelihood engines: Wald rows (logit-scale Wald for probabilities, exact for
# linear combinations; Monte Carlo from the MVN for ratio/product targets).
cjs_wald <- function(eng, tg, n_mvn) {
  th <- eng$theta; V <- eng$vcov; z <- stats::qnorm(0.975)
  mv <- NULL
  rows <- lapply(tg, function(g) {
    if (g$type == "coef") {
      se <- sqrt(V[g$idx, g$idx])
      return(data.frame(param = g$name, estimate = th[g$idx], sd = se, lower = th[g$idx] - z * se, upper = th[g$idx] + z * se))
    }
    if (g$type == "lc") {
      lp <- sum(th[g$idx] * g$x); se <- sqrt(max(drop(t(g$x) %*% V[g$idx, g$idx, drop = FALSE] %*% g$x), 0))
      e <- stats::plogis(lp)
      return(data.frame(param = g$name, estimate = e, sd = e * (1 - e) * se,
                        lower = stats::plogis(lp - z * se), upper = stats::plogis(lp + z * se)))
    }
    if (is.null(mv)) mv <<- cjs_mvn(th, V, n_mvn)
    v <- cjs_eval(mv, list(g))[, 1]
    e <- cjs_eval(matrix(th, 1, dimnames = list(NULL, names(th))), list(g))[1, 1]
    data.frame(param = g$name, estimate = e, sd = stats::sd(v), lower = stats::quantile(v, 0.025, names = FALSE),
               upper = stats::quantile(v, 0.975, names = FALSE))
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (!all(is.finite(V))) out$sd <- out$lower <- out$upper <- NA_real_
  out
}

cjs_fit <- function(data, settings) {
  settings <- merge_list(cjs_settings, settings)
  set.seed(settings$seed)
  engine <- match.arg(settings$engine, c("marked", "nimble"))
  des <- cjs_prepare(data, settings)
  eng <- switch(engine, marked = cjs_engine_marked(des, settings), nimble = cjs_engine_nimble(des, settings))
  tg <- cjs_targets(des)
  notes <- eng$notes
  gof <- NULL
  if (isTRUE(settings$gof)) {
    gof <- tryCatch(cjs_gof(data), error = function(e) { notes <<- c(notes, paste("R2ucare goodness-of-fit tests failed:", conditionMessage(e))); NULL })
    if (!is.null(gof)) { gof$model <- cjs_gof_model(gof, des); notes <- c(notes, cjs_gof_note(gof, des)) }
  }
  if (engine == "marked") {
    if (isTRUE(settings$chat_adjust) && !is.null(gof) && is.finite(gof$model$chat) && gof$model$chat > 1) {
      eng$vcov <- eng$vcov * gof$model$chat
      notes <- c(notes, sprintf("Standard errors inflated by sqrt(c-hat) = %.2f (quasi-likelihood adjustment for overdispersion).", sqrt(gof$model$chat)))
    }
    est <- cjs_wald(eng, tg, settings$n_mvn)
    dg <- list(converged = isTRUE(eng$converged), rhat_max = NA_real_, ess_min = NA_real_)
    # boundary check: a probability within 1e-4 of 0 or 1 is not estimable
    lc <- vapply(tg, function(g) g$type == "lc", TRUE)
    bnd <- est$param[lc][est$estimate[lc] > 0.999 | est$estimate[lc] < 0.001]
    cf <- vapply(tg, function(g) g$type == "coef", TRUE)
    wild <- est$param[cf][is.finite(est$sd[cf]) & (abs(est$estimate[cf]) > 10 | est$sd[cf] > 10)]
    if (length(bnd)) notes <- c(notes, sprintf("Estimates on the boundary (probability < 0.001 or > 0.999): %s. These are not estimable with these data; simplify the model (fewer terms, pool groups or occasions).", paste(bnd, collapse = ", ")))
    if (length(wild)) notes <- c(notes, sprintf("Coefficients with |estimate| > 10 or SE > 10 on the logit scale (not identifiable with these data): %s.", paste(wild, collapse = ", ")))
    if (length(bnd) || length(wild)) dg$converged <- FALSE
    draws <- if (all(is.finite(eng$vcov))) cjs_mvn(eng$theta, eng$vcov, settings$n_mvn) else NULL
  } else {
    pub <- lapply(eng$chains, function(TH) cjs_eval(TH, tg))
    s <- summarise_draws_eco(pub)
    est <- s$estimates
    dg <- list(converged = s$diagnostics$converged, rhat_max = s$diagnostics$rhat_max, ess_min = s$diagnostics$ess_min)
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      notes <- c(notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                dg$rhat_max, dg$ess_min))
    draws <- eng$draws
  }
  if (des$confounded_last)
    notes <- c(notes, "Survival and detection are both time-dependent: the last survival and last detection probability cannot be separated, so only their product (survival_x_detection_*) is reported (estimated by fixing recapture at the last occasion to 1).")
  if (!identical(des$formulas$phi, des$formulas$phi_user))
    notes <- c(notes, sprintf("Transients: the survival formula '%s' was fitted as '%s', so covariate and group effects (phi_beta_*) refer to residents (intervals after the first one following marking) and are not diluted by transients; phi_beta_tsm_first_* are first-interval deviations.",
                              des$formulas$phi_user, des$formulas$phi))
  if (des$random_time && engine == "nimble")
    notes <- c(notes, "Random time effect on survival: survival_* rows are survival in an average interval (random effect = 0, the median across intervals on the probability scale); survival_temporal_sd is the SD of the interval effects on the logit scale.")
  rownames(est) <- NULL
  dg$notes <- notes
  obj <- list(engine = eng$engine, des = des, targets = tg, theta = eng$theta, vcov = eng$vcov, draws = draws,
              chains = eng$chains, eps_chains = eng$eps_chains, loglik = eng$loglik, aic = eng$aic,
              gof = gof, settings = settings, time = eng$time)
  list(estimates = est, diagnostics = dg, object = obj)
}

# ---------------------------------------------------------------- goodness of fit
cjs_gof <- function(data) {
  eco_require("R2ucare")
  ch <- data$ch
  keep <- rowSums(ch) > 0
  grp <- if (!is.null(data$group)) as.character(data$ind[[data$group]]) else rep("all", nrow(ch))
  comps <- c("test3sr", "test3sm", "test2ct", "test2cl")
  lab <- c(test3sr = "3.SR (transience)", test3sm = "3.SM", test2ct = "2.CT (trap dependence)", test2cl = "2.CL")
  per <- list(); sign <- list()
  for (g in unique(grp)) {
    X <- ch[keep & grp == g, , drop = FALSE]
    if (nrow(X) < 5) next
    f <- rep(1, nrow(X))
    r <- list(test3sr = R2ucare::test3sr(X, f, verbose = FALSE), test3sm = R2ucare::test3sm(X, f, verbose = FALSE),
              test2ct = R2ucare::test2ct(X, f, verbose = FALSE), test2cl = R2ucare::test2cl(X, f, verbose = FALSE))
    for (k in comps) {
      v <- r[[k]][[1]]
      per[[length(per) + 1]] <- data.frame(group = g, test = lab[[k]], chi2 = as.numeric(v[["stat"]]), df = as.numeric(v[["df"]]))
      if (k %in% c("test3sr", "test2ct")) sign[[length(sign) + 1]] <- data.frame(group = g, test = k, z = as.numeric(v[["sign_test"]]))
    }
  }
  per <- do.call(rbind, per)
  tab <- stats::aggregate(cbind(chi2, df) ~ test, data = per, FUN = sum)
  tab <- tab[match(unname(lab), tab$test), ]; tab <- tab[!is.na(tab$test), ]
  tab$p_value <- ifelse(tab$df > 0, stats::pchisq(tab$chi2, tab$df, lower.tail = FALSE), NA_real_)
  tot <- data.frame(test = "Overall CJS (sum)", chi2 = sum(tab$chi2), df = sum(tab$df))
  tot$p_value <- if (tot$df > 0) stats::pchisq(tot$chi2, tot$df, lower.tail = FALSE) else NA_real_
  tab <- rbind(tab, tot); rownames(tab) <- NULL
  sg <- if (length(sign)) do.call(rbind, sign) else NULL
  zsum <- function(k) { if (is.null(sg)) return(NA_real_); z <- sg$z[sg$test == k & is.finite(sg$z)]; if (!length(z)) NA_real_ else sum(z) / sqrt(length(z)) }
  list(table = tab, chat = if (tot$df > 0) tot$chi2 / tot$df else NA_real_, p = tot$p_value,
       z_transience = zsum("test3sr"), z_trap = zsum("test2ct"), by_group = per,
       groups = unique(grp))
}

# c-hat for the FITTED model: components that the model already accounts for are removed
# from the overall test (3.SR when tsm is in the survival formula, 2.CT when td is in the
# detection formula; Pradel et al. 2005; Choquet et al. 2009 U-CARE guide).
cjs_gof_model <- function(gof, des) {
  t <- gof$table[gof$table$test != "Overall CJS (sum)", ]
  drop <- c(if ("tsm" %in% des$mphi$vars) "3.SR (transience)", if ("td" %in% des$mp$vars) "2.CT (trap dependence)")
  u <- t[!t$test %in% drop, ]
  chi2 <- sum(u$chi2); df <- sum(u$df)
  list(chi2 = chi2, df = df, chat = if (df > 0) chi2 / df else NA_real_,
       p = if (df > 0) stats::pchisq(chi2, df, lower.tail = FALSE) else NA_real_, excluded = drop,
       flagged = t$test[is.finite(t$p_value) & t$p_value < 0.05])
}

cjs_gof_note <- function(gof, des) {
  t <- gof$table
  g <- function(lbl) t[t$test == lbl, ]
  s <- sprintf("Goodness of fit (R2ucare, tests of the fully time-dependent CJS model%s): overall chi-square %.1f on %d df, p = %.3f, c-hat = %.2f.",
               if (length(gof$groups) > 1) ", summed over groups" else "", gof$table$chi2[nrow(t)], as.integer(t$df[nrow(t)]), gof$p, gof$chat)
  sr <- g("3.SR (transience)"); ct <- g("2.CT (trap dependence)")
  if (nrow(sr) && isTRUE(sr$p_value < 0.05)) {
    s <- c(s, sprintf("Test 3.SR detects transience (chi-square %.1f, df %d, p = %.3f%s): newly marked animals are seen again less often than previously marked ones. %s",
                      sr$chi2, as.integer(sr$df), sr$p_value,
                      if (is.finite(gof$z_transience)) sprintf(", signed z = %.2f", gof$z_transience) else "",
                      if ("tsm" %in% des$mphi$vars) "The survival formula includes tsm, which accounts for this."
                      else "Add 'tsm' to the survival formula (first interval after marking vs later); otherwise survival is biased low."))
  }
  if (nrow(ct) && isTRUE(ct$p_value < 0.05)) {
    s <- c(s, sprintf("Test 2.CT detects trap dependence (chi-square %.1f, df %d, p = %.3f%s; negative = trap-happy, positive = trap-shy). %s",
                      ct$chi2, as.integer(ct$df), ct$p_value,
                      if (is.finite(gof$z_trap)) sprintf(", signed z = %.2f", gof$z_trap) else "",
                      if ("td" %in% des$mp$vars) "The detection formula includes td, which accounts for this."
                      else "Add 'td' to the detection formula."))
  }
  gm <- gof$model
  if (length(gm$excluded))
    s <- c(s, sprintf("For the fitted model (which accounts for %s), the remaining components give chi-square %.1f on %d df, p = %.3f, c-hat = %.2f.",
                      paste(sub(" .*", "", gm$excluded), collapse = " and "), gm$chi2, as.integer(gm$df), gm$p, gm$chat))
  if (is.finite(gm$chat) && gm$chat > 3)
    s <- c(s, "c-hat > 3 indicates structural lack of fit, not just overdispersion: revise the model (transients, trap dependence, groups) before interpreting.")
  else if (is.finite(gm$chat) && gm$chat > 1.5)
    s <- c(s, "c-hat between 1.5 and 3: moderate overdispersion; consider settings chat_adjust: true (inflates standard errors by sqrt(c-hat)).")
  s
}
