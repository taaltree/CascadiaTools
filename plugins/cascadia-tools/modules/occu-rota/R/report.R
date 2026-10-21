# occu-rota: predictions, conditional-effects plots, interpretation, methods.

rota_fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
rota_ci <- function(e, l, u, d = 2) sprintf("%s (95%% CI %s to %s)", rota_fmt(e, d), rota_fmt(l, d), rota_fmt(u, d))
rota_cap <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
rota_lab <- function(x) gsub("_", " ", x)

# Biological-scale predictions over newdata rows (raw covariate units).
# Returns long data.frame: row, quantity, estimate, lower, upper.
rota_predict_state <- function(fit, newdata, prob = 0.95) {
  o <- fit$object; des <- o$des
  B <- o$draws[, seq_len(des$P), drop = FALSE]
  Xn <- rota_Xnew(des, newdata)
  Fl <- rota_F_from(des, B, Xn)
  Fp <- if (!is.null(o$point)) rota_F_from(des, matrix(o$point[seq_len(des$P)], 1), Xn) else NULL
  a <- (1 - prob) / 2
  out <- lapply(seq_len(nrow(Xn)), function(r) {
    D <- rota_derived(rota_psi(matrix(Fl[[r]], nrow(B)), des$A), des$Z)
    est <- if (!is.null(Fp)) as.vector(rota_derived(rota_psi(matrix(Fp[[r]], 1), des$A), des$Z)) else colMeans(D, na.rm = TRUE)
    data.frame(row = r, quantity = colnames(D), estimate = est,
               lower = apply(D, 2, stats::quantile, a, na.rm = TRUE),
               upper = apply(D, 2, stats::quantile, 1 - a, na.rm = TRUE), row.names = NULL)
  })
  do.call(rbind, out)
}

rota_predict_det <- function(fit, s, newdata, prob = 0.95) {
  o <- fit$object; des <- o$des
  si <- match(s, des$species); cols <- des$P + which(des$det_of == si)
  Wn <- rota_Wnew(des, s, newdata)
  Pd <- stats::plogis(Wn %*% t(o$draws[, cols, drop = FALSE]))
  a <- (1 - prob) / 2
  est <- if (!is.null(o$point)) as.vector(stats::plogis(Wn %*% o$point[cols])) else rowMeans(Pd)
  data.frame(estimate = est, lower = apply(Pd, 1, stats::quantile, a), upper = apply(Pd, 1, stats::quantile, 1 - a))
}

rota_state_vars <- function(des) unique(unlist(lapply(des$forms, function(f) all.vars(stats::as.formula(f)))))
rota_det_vars <- function(des) unique(unlist(lapply(des$dforms, function(f) all.vars(stats::as.formula(f)))))

rota_grid <- function(des, v, n = 60) {
  r <- range(des$raw[[v]], na.rm = TRUE)
  stats::setNames(data.frame(seq(r[1], r[2], length.out = n)), v)
}

rota_plots <- function(fit, data) {
  o <- fit$object; des <- o$des
  if (o$engine == "naive") return(list())
  sp <- des$species; pl <- list()
  # conditional occupancy at average conditions
  pr0 <- rota_predict_state(fit, data.frame())
  cond <- pr0[grepl("_given_", pr0$quantity), ]
  cond$target <- sub("^prob_(.*)_given_.*$", "\\1", cond$quantity)
  cond$other <- sub("^prob_.*_given_(.*)_(present|absent)$", "\\1", cond$quantity)
  cond$status <- ifelse(grepl("_present$", cond$quantity), "present", "absent")
  cond$lab <- paste0(rota_lab(cond$target), "\n(", rota_lab(cond$other), " ", cond$status, ")")
  marg <- pr0[grepl("^occupancy_prob_", pr0$quantity), ]
  marg$target <- sub("^occupancy_prob_", "", marg$quantity); marg$other <- "(any)"; marg$status <- "marginal"
  marg$lab <- paste0(rota_lab(marg$target), "\n(marginal)")
  cc <- rbind(cond, marg)
  cc$status <- factor(cc$status, levels = c("marginal", "absent", "present"))
  pl[["Occupancy of each species at average conditions: marginal, and conditional on another species"]] <-
    ggplot2::ggplot(cc, ggplot2::aes(x = estimate, y = lab, colour = status)) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.5) +
    ggplot2::scale_colour_manual(values = c(marginal = "grey35", absent = "#b5651d", present = "#2f6f4f")) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Probability the species occupies a site (95% interval)", y = NULL,
                  colour = "Other species") +
    ggplot2::theme_minimal(base_size = 12)
  # conditional effects along each numeric site covariate
  for (v in rota_state_vars(des)) {
    if (!is.null(des$levels[[v]])) {
      nd <- stats::setNames(data.frame(des$levels[[v]]), v)
      pr <- rota_predict_state(fit, nd); pr$level <- nd[[v]][pr$row]
      d <- pr[grepl("^occupancy_prob_", pr$quantity), ]; d$species <- rota_lab(sub("^occupancy_prob_", "", d$quantity))
      pl[[paste0("Occupancy by ", v)]] <- ggplot2::ggplot(d, ggplot2::aes(level, estimate, colour = species)) +
        ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), position = ggplot2::position_dodge(0.4)) +
        ggplot2::labs(x = v, y = "Marginal occupancy probability", colour = NULL) + ggplot2::theme_minimal(base_size = 12)
      next
    }
    g <- rota_grid(des, v)
    pr <- rota_predict_state(fit, g); pr$x <- g[[v]][pr$row]
    d <- pr[grepl("^occupancy_prob_", pr$quantity), ]; d$group <- rota_lab(sub("^occupancy_prob_", "", d$quantity))
    pl[[paste0("Marginal occupancy vs ", v)]] <- cond_effect_plot(
      d[, c("x", "estimate", "lower", "upper", "group")], v, "Occupancy probability (marginal)",
      observed_x = des$raw[[v]]) + ggplot2::coord_cartesian(ylim = c(0, 1))
    for (a in sp) for (b in sp) if (a != b) {
      q <- paste0("prob_", a, "_given_", b, c("_present", "_absent"))
      d <- pr[pr$quantity %in% q, ]
      d$group <- paste(rota_lab(b), ifelse(grepl("_present$", d$quantity), "present", "absent"))
      pl[[sprintf("%s occupancy vs %s, when %s is present or absent", rota_cap(rota_lab(a)), v, rota_lab(b))]] <-
        cond_effect_plot(d[, c("x", "estimate", "lower", "upper", "group")], v,
                         sprintf("P(%s present)", rota_lab(a)), observed_x = des$raw[[v]]) +
        ggplot2::coord_cartesian(ylim = c(0, 1))
    }
  }
  # detection
  for (v in rota_det_vars(des)) {
    if (!is.null(des$levels[[v]]) || is.null(des$raw[[v]])) next
    g <- rota_grid(des, v)
    d <- do.call(rbind, lapply(sp, function(s) {
      f <- rota_dforms_has(des, s, v); if (!f) return(NULL)
      cbind(x = g[[v]], rota_predict_det(fit, s, g), group = rota_lab(s))
    }))
    if (!is.null(d)) pl[[paste0("Detection probability per survey vs ", v)]] <-
      cond_effect_plot(d, v, "Detection probability per survey (where present)", observed_x = des$raw[[v]]) +
      ggplot2::coord_cartesian(ylim = c(0, 1))
  }
  pl
}

rota_dforms_has <- function(des, s, v) v %in% all.vars(stats::as.formula(des$dforms[[s]]))

rota_nice_unit <- function(sdv) {
  u <- 10^floor(log10(sdv)); c(u, 2 * u, 5 * u)[which.min(abs(c(u, 2 * u, 5 * u) - sdv))]
}

rota_interpret <- function(fit, data) {
  o <- fit$object; des <- o$des; est <- fit$estimates; sp <- des$species
  g <- function(p) est[est$param == p, , drop = FALSE]
  out <- character()
  if (o$engine == "naive")
    return("WARNING: this is the naive control analysis that ignores imperfect detection. Its estimates are biased and must not be reported.")
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge cleanly or some parameters are at a boundary (see Diagnostics). Do not interpret the numbers below until this is resolved.")
  y0 <- des$y; y0[is.na(y0)] <- 0
  det <- apply(y0, c(1, 3), max)
  out <- c(out, sprintf("Data: %d sites, median %d surveys per site. Sites with at least one detection: %s. %s",
                        des$M, as.integer(stats::median(rowSums(des$w))),
                        paste(sprintf("%s %d (naive occupancy %.2f)", rota_lab(sp), colSums(det), colMeans(det)), collapse = "; "),
                        paste(vapply(utils::combn(des$S, 2, simplify = FALSE), function(cc)
                          sprintf("Both %s and %s were detected at %d sites.", rota_lab(sp[cc[1]]), rota_lab(sp[cc[2]]),
                                  sum(det[, cc[1]] * det[, cc[2]])), ""), collapse = " ")))
  avg <- if (length(rota_state_vars(des))) " with average covariate values" else ""
  hv <- unique(c(rota_state_vars(des), rota_det_vars(des)))
  if (length(hv)) out <- c(out, sprintf("'Average conditions' / 'average covariate values' below means a site and survey with %s. Every curve holds the covariates it does not vary at these values.",
                                        rota_hold_text(des, hv)))
  avg_lead <- if (nzchar(avg)) "At average covariate values, " else ""
  m <- vapply(sp, function(s) { r <- g(paste0("occupancy_prob_", s)); sprintf("%s %s", rota_lab(s), rota_ci(r$estimate, r$lower, r$upper)) }, "")
  out <- c(out, sprintf("Allowing for missed detections, the estimated probability that each species occupies a site%s is: %s.",
                        avg, paste(m, collapse = "; ")))
  # pairwise conditional occupancy and interaction
  for (cc in utils::combn(des$S, 2, simplify = FALSE)) {
    a <- sp[cc[1]]; b <- sp[cc[2]]
    ip <- g(paste0("interaction_", a, "_", b))
    if (!nrow(ip)) { out <- c(out, sprintf("The %s-%s interaction was fixed at zero (formula \"0\"), so the model assumes they occur independently.", rota_lab(a), rota_lab(b))); next }
    ga <- g(paste0("prob_", b, "_given_", a, "_present")); gn <- g(paste0("prob_", b, "_given_", a, "_absent"))
    ha <- g(paste0("prob_", a, "_given_", b, "_present")); hn <- g(paste0("prob_", a, "_given_", b, "_absent"))
    or <- exp(c(ip$estimate, ip$lower, ip$upper))
    dir <- if (ip$lower > 0) sprintf("more often than expected if they occurred independently (positive co-occurrence)")
           else if (ip$upper < 0) sprintf("less often than expected if they occurred independently (negative co-occurrence; this pattern alone does not show that one species avoids the other)")
           else "about as often as expected if they occurred independently: the data cannot distinguish a positive from a negative association"
    third <- if (des$S == 3) sprintf(" (for three species, this odds ratio applies where %s is absent; the conditional probabilities average over it)", rota_lab(setdiff(sp, c(a, b)))) else ""
    out <- c(out, sprintf(
      "%s%s and %s occur together %s. Where %s is present, %s occupies a site with probability %s; where %s is absent, %s. Equivalently, P(%s | %s present) = %s vs P(%s | %s absent) = %s. The co-occurrence odds ratio is %s%s; a value of 1 would mean no association.",
      avg_lead, if (nzchar(avg_lead)) rota_lab(a) else rota_cap(rota_lab(a)), rota_lab(b), dir,
      rota_lab(a), rota_lab(b), rota_ci(ga$estimate, ga$lower, ga$upper), rota_lab(a), rota_ci(gn$estimate, gn$lower, gn$upper),
      rota_lab(a), rota_lab(b), rota_ci(ha$estimate, ha$lower, ha$upper), rota_lab(a), rota_lab(b), rota_ci(hn$estimate, hn$lower, hn$upper),
      rota_ci(or[1], or[2], or[3]), third))
  }
  # covariate effects on occupancy and on the interaction
  for (v in rota_state_vars(des)) {
    if (!is.null(des$levels[[v]])) {
      out <- c(out, sprintf("'%s' is a categorical covariate; see the figure 'Occupancy by %s' for each level (reference level: %s).", v, v, des$levels[[v]][1]))
      next
    }
    qs <- stats::quantile(des$raw[[v]], c(0.1, 0.9))
    pr <- rota_predict_state(fit, stats::setNames(data.frame(qs), v))
    gp <- function(q, r) pr[pr$quantity == q & pr$row == r, ]
    txt <- vapply(sp, function(s) {
      q <- paste0("occupancy_prob_", s); lo <- gp(q, 1); hi <- gp(q, 2)
      sprintf("%s from %s to %s", rota_lab(s), rota_ci(lo$estimate, lo$lower, lo$upper), rota_ci(hi$estimate, hi$lower, hi$upper))
    }, "")
    out <- c(out, sprintf("Across the observed range of %s (10th percentile %s to 90th percentile %s), marginal occupancy changes: %s. This is an association, not necessarily a causal effect.",
                          v, signif(qs[1], 3), signif(qs[2], 3), paste(txt, collapse = "; ")))
    for (cc in utils::combn(des$S, 2, simplify = FALSE)) {
      a <- sp[cc[1]]; b <- sp[cc[2]]
      xs <- des$xscale[[paste0(a, ":", b)]]
      if (is.null(xs)) next
      for (cn in colnames(xs)[-1]) {
        cv <- tryCatch(all.vars(str2lang(cn)), error = function(e) character())
        if (length(cv) != 1 || cv != v) next
        r <- g(paste0("interaction_", a, "_", b, "_beta_", rota_public_col(cn)))
        if (!nrow(r)) next
        sdv <- xs["sd", cn]; u <- rota_nice_unit(sdv)
        lo <- gp(paste0("prob_", b, "_given_", a, "_present"), 1); hi <- gp(paste0("prob_", b, "_given_", a, "_present"), 2)
        lo0 <- gp(paste0("prob_", b, "_given_", a, "_absent"), 1); hi0 <- gp(paste0("prob_", b, "_given_", a, "_absent"), 2)
        out <- c(out, sprintf(
          "The %s-%s co-occurrence odds ratio is multiplied by %s per 1 SD (%s) increase in %s%s%s. At low %s (10th percentile), P(%s | %s present) = %s vs %s absent %s; at high %s (90th percentile), %s vs %s.",
          rota_lab(a), rota_lab(b), rota_ci(exp(r$estimate), exp(r$lower), exp(r$upper)), signif(sdv, 3), cn,
          if (cn == v && abs(u / sdv - 1) > 0.05) sprintf(", i.e. by %s per %s units", rota_fmt(exp(r$estimate * u / sdv)), signif(u, 3))
          else if (grepl(paste0("^log\\(\\s*", v, "\\s*\\)$"), cn)) sprintf(", i.e. by %s (95%% CI %s to %s) per doubling of %s",
            rota_fmt(exp(r$estimate * log(2) / sdv)), rota_fmt(exp(r$lower * log(2) / sdv)), rota_fmt(exp(r$upper * log(2) / sdv)), v) else "",
          if (r$lower > 0 || r$upper < 0) "" else " (the interval includes 1, so the data do not show that the association changes along this gradient)",
          v, rota_lab(b), rota_lab(a), rota_ci(lo$estimate, lo$lower, lo$upper), rota_lab(a), rota_ci(lo0$estimate, lo0$lower, lo0$upper),
          v, rota_ci(hi$estimate, hi$lower, hi$upper), rota_ci(hi0$estimate, hi0$lower, hi0$upper)))
      }
    }
  }
  # detection
  dt <- vapply(sp, function(s) { r <- g(paste0("detection_prob_", s)); sprintf("%s %s", rota_lab(s), rota_ci(r$estimate, r$lower, r$upper)) }, "")
  out <- c(out, sprintf("Per-survey detection probability where the species is present (average survey conditions): %s. With %d surveys, the chance of detecting a species at least once at an occupied site is %s.",
                        paste(dt, collapse = "; "), as.integer(stats::median(rowSums(des$w))),
                        paste(vapply(sp, function(s) { p <- g(paste0("detection_prob_", s))$estimate
                          sprintf("%s %.2f", rota_lab(s), 1 - (1 - p)^stats::median(rowSums(des$w))) }, ""), collapse = ", ")))
  db <- est[grepl("^p_.*_beta_", est$param), ]
  if (nrow(db)) {
    info <- do.call(rbind, lapply(seq_along(des$wscale), function(s) {
      w <- des$wscale[[s]][, -1, drop = FALSE]
      if (!ncol(w)) return(NULL)
      data.frame(param = paste0("p_", des$species[s], "_beta_", rota_public_col(colnames(w))),
                 lab = paste0(rota_lab(des$species[s]), ", ", colnames(w)),
                 dummy = w["mean", ] == 0 & w["sd", ] == 1)
    }))
    k <- match(db$param, info$param)
    unit <- ifelse(info$dummy[k] %in% TRUE, "vs the reference level", sprintf("per 1 SD"))
    ttx <- sprintf("%s: odds ratio %s %s", info$lab[k], unit, rota_ci(exp(db$estimate), exp(db$lower), exp(db$upper)))
    out <- c(out, paste0("Detection covariates (these describe how easy a species is to detect, not where it lives): ", paste(ttx, collapse = "; "), "."))
  }
  out <- c(out, "Caution: co-occurrence is not proof of a species interaction. A positive or negative association can arise from shared or opposite responses to habitat that the model does not include, from different spatial scales of use, or from one species affecting the other's detectability rather than its occupancy (Blanchet et al. 2020). Interpret the interaction as 'occur together more/less than expected given the covariates in the model'.")
  gof <- fit$diagnostics$gof
  if (!is.null(gof)) out <- c(out, if (gof$p_value < 0.05)
    sprintf("Goodness of fit: the model does NOT adequately reproduce the detection data (p = %.3f, c-hat %.2f). Intervals are probably too narrow; consider additional covariates or a different model before relying on these results.", gof$p_value, gof$c_hat)
    else sprintf("Goodness of fit: no evidence of lack of fit (%s, p = %.2f, c-hat %.2f).", gof$method, gof$p_value, gof$c_hat))
  out
}

rota_methods_text <- function(fit, data, settings) {
  o <- fit$object; des <- o$des; st <- o$settings
  fs <- paste(sprintf("%s: %s", names(des$forms), unlist(des$forms)), collapse = "; ")
  dfs <- paste(sprintf("%s: %s", names(des$dforms), unlist(des$dforms)), collapse = "; ")
  zero <- setdiff(names(des$terms_all), names(des$terms))
  eng <- switch(o$engine,
    unmarked = sprintf("Parameters were estimated by maximum likelihood with occuMulti in the R package unmarked (Fiske & Chandler 2011; Kellner et al. 2023)%s. Uncertainty in occupancy, conditional occupancy and co-occurrence probabilities was propagated by drawing %d parameter vectors from the multivariate normal approximation to the sampling distribution of the estimates; 95%% intervals are the 2.5%% and 97.5%% quantiles.",
                       if (isTRUE(st$penalty > 0)) sprintf(" with the penalized likelihood of Clipp et al. (2021), penalty %.2f, standard errors from %d bootstrap samples", st$penalty, st$penalty_boot) else "",
                       st$n_draws),
    nimble = sprintf("The model was fitted in a Bayesian framework with NIMBLE (de Valpine et al. 2017), using the categorical latent-state formulation: each site's occupancy state is one of the %d species combinations, z_i ~ Categorical(psi_i). Priors: Normal(0, SD %g) on every natural-parameter and logit-detection coefficient. %d chains of %d iterations (%d burn-in); convergence was assessed with R-hat < 1.01 and bulk effective sample size >= 400.",
                     des$C, st$prior_sd, st$nchains, st$niter, st$nburnin),
    naive = "NAIVE CONTROL: the multivariate Bernoulli model was fitted to 'detected at least once' as if detection were perfect.")
  paste0(
    sprintf("We modelled the co-occurrence of %s at %d sites with the multispecies occupancy model of Rota et al. (2016). ", { sl <- rota_lab(des$species); if (length(sl) > 1) paste(paste(sl[-length(sl)], collapse = ", "), sl[length(sl)], sep = " and ") else sl }, des$M),
    sprintf("The latent occupancy state of the %d species at each site was a multivariate Bernoulli random variable whose %d state probabilities are a multinomial-logit function of first-order natural parameters (one per species) and %s natural parameters (log odds ratios of co-occurrence)%s. ",
            des$S, des$C, if (des$max_order == 3) "second- and third-order" else "second-order",
            if (length(zero)) paste0("; ", paste(zero, collapse = ", "), " fixed at 0") else ""),
    sprintf("Natural parameters were linear functions of site covariates (%s). Detections on each survey were independent Bernoulli outcomes given presence, with species-specific detection probability on the logit scale (%s). ", fs, dfs),
    if (isTRUE(des$standardize)) "Numeric covariates were standardised (mean 0, SD 1) after any transformation, so intercepts refer to average conditions and coefficients are per 1 SD. " else "",
    eng,
    " Marginal occupancy, co-occurrence probability and the probability that each species occupies a site conditional on another species being present or absent were derived at average covariate values (the covariate value at which each standardised model term is at its mean, e.g. the geometric mean for log-transformed covariates; factors at their reference level) and along each covariate gradient. ",
    "Goodness of fit was assessed with a chi-square statistic on the frequency of sites with 0, 1, 2, ... detections of each species plus a chi-square statistic on the number of sites where each pair was detected together, by parametric bootstrap (maximum likelihood) or posterior predictive check (Bayesian). ",
    "The estimator was validated by simulation (ecotools module occu-rota; see VALIDATION.md for recovery, coverage and false-positive rates).")
}
