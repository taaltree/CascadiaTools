# scr: design notes, interpretation, figures, methods text, references.

scr_references <- c(
  "Efford M (2004) Density estimation in live-trapping studies. Oikos 106:598-610. doi:10.1111/j.0030-1299.2004.13043.x",
  "Borchers DL, Efford MG (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics 64:377-385. doi:10.1111/j.1541-0420.2007.00927.x",
  "Royle JA, Young KV (2008) A hierarchical model for spatial capture-recapture data. Ecology 89:2281-2289. doi:10.1890/07-0601.1",
  "Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press, Waltham MA.",
  "Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152(Suppl 2):S521-S537. doi:10.1007/s10336-010-0619-4",
  "Sollmann R, Gardner B, Belant JL (2012) How does spatial study design influence density estimates from spatial capture-recapture models? PLoS ONE 7:e34575. doi:10.1371/journal.pone.0034575",
  "Sun CC, Fuller AK, Royle JA (2014) Trap configuration and spacing influences parameter estimates in spatial capture-recapture models. PLoS ONE 9:e88025. doi:10.1371/journal.pone.0088025",
  "Milleret C, Dupont P, Bonenfant C, Broseth H, Flagstad O, Sutherland C, Bischof R (2019) A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363. doi:10.1002/ece3.4751",
  "Efford MG (2025) secr: Spatially explicit capture-recapture models. R package version 5.x. doi:10.32614/CRAN.package.secr",
  "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487",
  "Efford MG, Borchers DL, Mowat G (2013) Varying effort in capture-recapture studies. Methods in Ecology and Evolution 4:629-636. doi:10.1111/2041-210X.12049",
  "Gardner B, Royle JA, Wegan MT, Rainbolt RE, Curtis PD (2010) Estimating black bear density using DNA data from hair snares. Journal of Wildlife Management 74:318-325. doi:10.2193/2009-101",
  "Gardner B, Royle JA, Wegan MT (2009) Hierarchical models for estimating density from DNA mark-recapture studies. Ecology 90:1106-1115.",
  "Pledger S (2000) Unified maximum likelihood estimates for closed capture-recapture models using mixtures. Biometrics 56:434-442. doi:10.1111/j.0006-341X.2000.00434.x",
  "Royle JA, Gardner B (2011) Hierarchical spatial capture-recapture models for estimating density from trapping arrays. In: O'Connell AF, Nichols JD, Karanth KU (eds) Camera Traps in Animal Ecology. Springer, Tokyo, pp 163-190.",
  "Lukacs PM, Burnham KP (2005) Review of capture-recapture methods applicable to noninvasive genetic sampling. Molecular Ecology 14:3909-3919. doi:10.1111/j.1365-294X.2005.02717.x",
  "Paetkau D (2003) An empirical exploration of data quality in DNA-based population inventories. Molecular Ecology 12:1375-1387.",
  "Taberlet P, Waits LP, Luikart G (1999) Noninvasive genetic sampling: look before you leap. Trends in Ecology & Evolution 14:323-327.",
  "Link WA, Yoshizaki J, Bailey LL, Pollock KH (2010) Uncovering a latent multinomial: analysis of mark-recapture data with misidentification. Biometrics 66:178-185. doi:10.1111/j.1541-0420.2009.01244.x",
  "Phan D, Pradhan N, Jankowiak M (2019) Composable effects for flexible and accelerated probabilistic programming in NumPyro. arXiv:1912.11554")

# ---------------------------------------------------------------- draws
# Common parameter draws for every engine: psi, det0, sigma, beta*, phi
# (internal units) with attribute M. For secr, pseudo-draws from the
# multivariate normal of the link-scale MLEs, expressed as psi = E[N] / M with
# large M.
scr_pdraws <- function(fit, ndraw = 2000) {
  ob <- fit$object
  if (!is.null(ob$chains)) {
    d <- scr_pool(ob$chains, ndraw); attr(d, "M") <- ob$prep$M_used; return(d)
  }
  prep <- ob$prep
  B <- MASS_mvrnorm(ndraw, ob$secr_b, ob$secr_V); colnames(B) <- names(ob$secr_b)
  out <- t(apply(B, 1, function(bb) scr_secr_to_th(prep, stats::setNames(bb, names(ob$secr_b)), ob$secr_Mbig)))
  attr(out, "M") <- ob$secr_Mbig
  out
}

# Expected density (reported units) at standardised density-covariate vector xs for each draw
scr_density_at <- function(prep, pd, xs) {
  M <- attr(pd, "M"); PD <- ncol(prep$XD)
  vapply(seq_len(nrow(pd)), function(r) {
    b <- if (PD) pd[r, paste0("beta", seq_len(PD))] else numeric(0)
    lw <- if (PD) as.vector(prep$XD %*% b) else rep(0, prep$G)
    pd[r, "psi"] * M * exp(sum(b * xs) - lse(lw)) / prep$pix^2 * scr_area_conv(prep)
  }, 1)
}

scr_q <- function(x, p = c(0.025, 0.975)) stats::quantile(x, p, na.rm = TRUE, names = FALSE)
scr_lab <- function(v, data, opts = NULL) {
  d <- if (!is.null(data$covariate_labels[[v]])) data$covariate_labels[[v]] else v
  fig_label(v, d, opts)
}
# significant digits, keeping trailing zeros ("0.10") but never a bare trailing dot ("1325.")
scr_fmt <- function(x, d = 2) sub("\\.$", "", trimws(formatC(x, digits = d, format = "fg", flag = "#")))
scr_unit_len <- function(prep) if (is.na(prep$units)) "map units" else prep$units

# ---------------------------------------------------------------- design notes
scr_design_notes <- function(prep, est) {
  g <- function(p) est$estimate[est$param == p]
  sig <- est$estimate[est$param %in% c("sigma", "sigma_female", "sigma_male")]; out <- character()
  if (!length(sig) || !all(is.finite(sig))) return(out)
  sig <- max(sig)
  ratio <- prep$trap_spacing / sig
  out <- c(out, sprintf("Trap spacing (%s %s) is %.1f x sigma. Spacing of about 1-2 sigma or less is recommended (Sun et al. 2014); wider spacing yields few spatial recaptures and imprecise or biased sigma and density.",
                        format(signif(prep$trap_spacing, 3)), scr_unit_len(prep), ratio))
  if (ratio > 2) out <- c(out, "WARNING: traps are more than 2 sigma apart; treat sigma and density with caution (see VALIDATION.md, harder design).")
  if (!is.null(prep$buffer) && is.finite(prep$buffer)) {
    br <- prep$buffer / sig
    out <- c(out, sprintf("State-space buffer is %.1f x sigma beyond the outermost traps.", br))
    if (br < 3) out <- c(out, "WARNING: the buffer is less than 3 sigma: animals living outside the state space can be detected and are squeezed into it, inflating density. Use a buffer of at least 3-4 sigma.")
  }
  if (prep$spacing / prep$cscale > sig / prep$cscale) out <- c(out, "WARNING: habitat-mask pixels are larger than sigma; use a finer mask (spacing ~ sigma/2).")
  out <- c(out, sprintf("%d of %d detected individuals were detected at more than one trap (spatial recaptures).", prep$n_spatial_recap, prep$n))
  if (prep$n_spatial_recap < 10) out <- c(out, sprintf("WARNING: only %d spatial recaptures (< 10): sigma, p0 and density are weakly identified. In simulations with traps 3 sigma apart (median 3 spatial recaptures) density was biased about +18%% and 14%% of fits did not converge (VALIDATION.md section 2).", prep$n_spatial_recap))
  dr <- est[est$param == scr_density_name(prep$units), , drop = FALSE]
  if ("rhat" %in% names(est) && nrow(dr) == 1 && is.finite(dr$sd) && dr$estimate > 0) {
    cv <- dr$sd / dr$estimate
    if (cv > 0.2) out <- c(out, sprintf("NOTE: density is imprecise (posterior SD = %.0f%% of the posterior mean). The posterior mean of a right-skewed posterior sits above the most likely value by roughly CV^2 (here about %.0f%%); in simulations at posterior CVs of 0.2-0.3 posterior-mean density was 4-10%% above secr's maximum-likelihood estimate and 9-20%% above the truth (VALIDATION.md, Independent verification; +19%% with ~46 detected animals in a sex-specific Poisson hair-snare design). Report the 95%% interval, and the posterior median or the secr estimate alongside the mean.",
                                         100 * cv, 100 * cv^2))
  }
  if (prep$obs != "poisson" && prep$behav == "trap") out <- c(out, "NOTE: with 0/1 detections, baseline detection p0 and the trap-specific response trade off: in simulations (VALIDATION.md 10f, 12c) posterior-mean p0 was ~16% too high and the response biased toward 'trap-shy' (-0.38 when the true response was 0), while density, abundance and sigma were unbiased. secr (maximum likelihood) was unbiased for both on the same design, so report its p0 and response (compare_engines: [secr]) and compare with the model without the response.")
  if (prep$n < 25) out <- c(out, sprintf("WARNING: only %d individuals detected (< 25). In simulations with ~18 detected animals (traps 1.5 sigma apart) density was biased about +12%% and p0 about +21%% (VALIDATION.md section 2b): treat density as an upper-leaning estimate and report its interval.", prep$n))
  out
}

# ---------------------------------------------------------------- interpret
scr_effect_sentence <- function(r, sc, lab, kind, pois) {
  scale_word <- if (pois) "the expected number of samples" else "the odds of detection"
  what <- if (kind == "sigma") "the spatial scale sigma" else scale_word
  if (isTRUE(sc$binary))
    sprintf("'%s' (1 vs 0) multiplied %s by %s (95%% interval %s to %s).", lab, what, scr_fmt(exp(r$estimate)), scr_fmt(exp(r$lower)), scr_fmt(exp(r$upper)))
  else sprintf("each 1 SD (%s) increase in '%s' multiplied %s by %s (95%% interval %s to %s).", scr_fmt(sc$scale), lab, what,
               scr_fmt(exp(r$estimate)), scr_fmt(exp(r$lower)), scr_fmt(exp(r$upper)))
}

scr_interpret <- function(fit, data) {
  est <- fit$estimates; prep <- fit$object$prep
  row <- function(p) est[est$param == p, ]
  pois <- prep$obs == "poisson"
  dn <- scr_density_name(prep$units); al <- scr_area_label(prep$units); detn <- scr_det_name(prep)
  D <- row(dn); N <- row("abundance")
  ci <- function(r, d = 2) sprintf("%s (95%% %s %s to %s)", scr_fmt(r$estimate, d), if (fit$object$engine == "secr") "CI" else "CrI",
                                   scr_fmt(r$lower, d), scr_fmt(r$upper, d))
  ndet <- if (pois) sum(prep$y) else nrow(prep$captures)
  A_rep <- prep$area * scr_hr_conv(prep)
  obs_txt <- switch(prep$obs,
    binomial = "0/1 per animal, detector and occasion, summed over the occasions each detector operated (binomial model)",
    bernoulli = "0/1 per animal, detector and occasion (Bernoulli model, occasion by occasion)",
    poisson = "the number of samples per animal, detector and occasion (Poisson model: one animal can leave several hair or scat samples at one detector)")
  out <- c(
    sprintf("Data: %d individuals were detected %d times (%s) at %d detectors over %d occasions; %d individuals were detected at more than one detector, which is what lets the model estimate how far animals range. Detections were modelled as %s.",
            prep$n, ndet, if (pois) "genetic samples" else "individual x detector x occasion", prep$J, prep$K, prep$n_spatial_recap, obs_txt),
    sprintf("Density was estimated at %s animals per %s%s. Over the %s %s of habitat in the state space (the detector array plus a buffer), that is %s animals, of whom %d (%.0f%%) were detected.",
            ci(D, 2), al, if (dn == "density_per_km2") sprintf(", i.e. %s per 100 km2", scr_fmt(100 * D$estimate, 3)) else "",
            scr_fmt(A_rep, 3), al, ci(N, 3), prep$n, 100 * prep$n / N$estimate))
  det_txt <- function(r, who = "an animal's") if (pois)
    sprintf("A detector placed at %s activity centre would collect on average %s samples of it per occasion with full effort (the baseline encounter rate lambda0)", who, ci(r, 2))
    else sprintf("A detector placed at %s activity centre would detect it on a given occasion with probability %s", who, ci(r, 2))
  ref <- if (length(c(prep$pcov, prep$ocov))) ", at the reference level of the detection covariates (0 for 0/1 covariates, the mean for continuous ones)" else ""
  if (prep$S == 1) {
    sg <- row("sigma"); hr <- row("home_range_area_95")
    out <- c(out,
      sprintf("The spatial scale of detection, sigma, was %s %s. If space use is a circular bivariate normal around each activity centre, 95%% of an animal's activity falls within %s %s of its centre, a home-range area of about %s %s. Real home ranges are rarely circular, so treat this as an index of ranging scale.",
              ci(sg, 3), scr_unit_len(prep), scr_fmt(sqrt(stats::qchisq(0.95, 2)) * sg$estimate, 3), scr_unit_len(prep), ci(hr, 2), al),
      paste0(det_txt(row(detn)), ref, "; detection falls off with distance (see the detection-function figure)."))
  } else {
    pf <- row("proportion_female")
    sf <- row("sigma_female"); sm <- row("sigma_male"); hf <- row("home_range_area_95_female"); hm <- row("home_range_area_95_male")
    out <- c(out,
      sprintf("Sex: the proportion of females in the population was estimated at %s; animals of unknown sex (%d of %d detected) were included with their sex treated as unknown.",
              ci(pf, 2), sum(prep$sexi < 0), prep$n),
      sprintf("Females: sigma %s %s, 95%% home range %s %s. Males: sigma %s %s, 95%% home range %s %s.%s",
              ci(sf, 3), scr_unit_len(prep), ci(hf, 2), al, ci(sm, 3), scr_unit_len(prep), ci(hm, 2), al,
              if ("sigma" %in% prep$sexeff) { r <- row("sigma_beta_male"); sprintf(" Males ranged %s times as far as females (95%% interval %s to %s).", scr_fmt(exp(r$estimate)), scr_fmt(exp(r$lower)), scr_fmt(exp(r$upper))) } else " (sigma was assumed equal for the sexes)"),
      paste0(det_txt(row(paste0(detn, "_female")), "a female's"), ref, "; for males: ", ci(row(paste0(detn, "_male")), 2), "."))
  }
  pd <- tryCatch(scr_pdraws(fit, 1000), error = function(e) NULL)
  for (k in seq_along(c(prep$pcov, prep$ocov))) {
    v <- c(prep$pcov, prep$ocov)[k]; r <- row(paste0("detection_beta_", v))
    sc <- if (v %in% prep$pcov) prep$pscale[[v]] else prep$oscale[[v]]
    lab <- scr_lab(v, data)
    lvl <- if (v %in% prep$pcov) "detector-level" else "occasion-level"
    s1 <- paste0("Detection (", lvl, "): ", scr_effect_sentence(r, sc, lab, "det", pois))
    if (isTRUE(sc$binary) && !is.null(pd) && prep$S == 1) {
      bi <- paste0("beta", prep$bidx[3] + k)
      p1 <- if (pois) pd[, "det0"] * exp(pd[, bi]) else stats::plogis(stats::qlogis(pd[, "det0"]) + pd[, bi])
      s1 <- paste(s1, sprintf("Baseline %s was %.2f (%.2f to %.2f) where %s = 1 versus %s where it = 0.", if (pois) "encounter rate" else "detection",
                              mean(p1), scr_q(p1)[1], scr_q(p1)[2], lab, scr_fmt(row(detn)$estimate)))
    }
    out <- c(out, paste(s1, "This is a detection effect (how easily animals are recorded), not a difference in where animals live."))
  }
  for (v in prep$scov) {
    r <- row(paste0("sigma_beta_", v))
    out <- c(out, paste0("Spatial scale: ", scr_effect_sentence(r, prep$sscale[[v]], scr_lab(v, data), "sigma", pois),
                         " A larger sigma at a detector means it draws in animals from farther away (e.g. a scent lure); it describes detection, not home-range size."))
  }
  if (prep$S == 2 && "p0" %in% prep$sexeff) {
    r <- row("detection_beta_male")
    out <- c(out, sprintf("Sex and detection: males had %s times the %s of females at the same distance (95%% interval %s to %s).",
                          scr_fmt(exp(r$estimate)), if (pois) "expected number of samples" else "odds of detection", scr_fmt(exp(r$lower)), scr_fmt(exp(r$upper))))
  }
  if (prep$behav == "trap") {
    r <- row("detection_beta_trap_response"); e <- exp(c(r$estimate, r$lower, r$upper))
    out <- c(out, sprintf("Behavioural response: after an animal was first detected at a detector, its %s at that detector %s multiplied by %s (95%% interval %s to %s): %s. In simulations with 0/1 data, ignoring a real trap-happy response mainly biased baseline detection upward (+43%%) and sigma downward (-6%%), with little effect on density (VALIDATION.md 12b); compare with the model without the response.",
                          if (pois) "expected number of samples" else "odds of detection", if (pois) "was" else "were", scr_fmt(e[1]), scr_fmt(e[2]), scr_fmt(e[3]),
                          if (e[2] > 1) "animals were 'trap-happy' (e.g. returning to a baited site)" else if (e[3] < 1) "animals were 'trap-shy'" else "the interval includes no response (1)"))
  }
  for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]; r <- row(paste0("density_beta_", v)); sc <- prep$dscale[[v]]; lab <- scr_lab(v, data)
    rr <- exp(c(r$estimate, r$lower, r$upper))
    s1 <- sprintf("Density and %s: each 1 SD (%s units) higher %s is associated with %s%% %s density (rate ratio %s, 95%% interval %s to %s).",
                  lab, scr_fmt(sc$scale), lab, scr_fmt(abs(100 * (rr[1] - 1)), 2), if (rr[1] >= 1) "higher" else "lower",
                  scr_fmt(rr[1]), scr_fmt(rr[2]), scr_fmt(rr[3]))
    if (!is.null(pd)) {
      xr <- stats::quantile(prep$mask[[v]], c(0.1, 0.9))
      xs <- function(x) { z <- rep(0, length(prep$dcov)); z[k] <- (x - sc$center) / sc$scale; z }
      lo <- scr_density_at(prep, pd, xs(xr[1])); hi <- scr_density_at(prep, pd, xs(xr[2]))
      pb <- mean(pd[, paste0("beta", k)] > 0)
      s1 <- paste(s1, sprintf("Expected density is %s per %s (%s to %s) where %s = %s (10th percentile of the habitat) and %s (%s to %s) where it = %s (90th percentile)%s.",
                              scr_fmt(mean(lo)), al, scr_fmt(scr_q(lo)[1]), scr_fmt(scr_q(lo)[2]), lab, scr_fmt(xr[1], 3),
                              scr_fmt(mean(hi)), scr_fmt(scr_q(hi)[1]), scr_fmt(scr_q(hi)[2]), scr_fmt(xr[2], 3),
                              if (fit$object$engine == "secr") "" else sprintf(", other density covariates at their mean. Posterior probability that the association is positive: %.2f", pb)))
    }
    s1 <- paste(s1, "This is a correlation between density and habitat across the state space, not proof that the habitat causes higher density.")
    out <- c(out, s1)
  }
  sgm <- max(est$estimate[est$param %in% c("sigma", "sigma_female", "sigma_male")])
  ratio <- prep$trap_spacing / sgm
  out <- c(out, sprintf("Design: detectors were %.1f sigma apart (recommended <= 2 sigma; Sun et al. 2014) and the state space extends %.1f sigma beyond the detectors (recommended >= 3-4 sigma)%s.",
                        ratio, if (!is.null(prep$buffer)) prep$buffer / sgm else NA, if (prep$S == 2) ", using the larger (male) sigma" else ""))
  bc <- fit$object$buffer_check
  if (!is.null(bc) && nrow(bc) == 2) {
    ch <- bc$estimate[2] / bc$estimate[1] - 1
    rr <- if (!is.null(bc$region_buffer)) sprintf(" within %s %s of the detectors (the region both fits share)", format(signif(bc$region_buffer[1], 3)), scr_unit_len(prep)) else ""
    out <- c(out, sprintf("Buffer check: refitting with a %s state space (%.1f instead of %.1f sigma) gave density%s of %s per %s (%s to %s), %+.1f%% from the main fit (%s): %s",
                          if (bc$buffer[2] > bc$buffer[1]) "larger" else "smaller", bc$buffer_sigma[2], bc$buffer_sigma[1], rr,
                          scr_fmt(bc$estimate[2]), al, scr_fmt(bc$lower[2]), scr_fmt(bc$upper[2]), 100 * ch, scr_fmt(bc$estimate[1]),
                          if (abs(ch) < 0.05 && abs((bc$estimate[2] - bc$estimate[1]) / bc$sd[1]) < 0.3) "the buffer is large enough." else "density depends on the buffer, so the smaller buffer is inadequate."))
  }
  if (!is.null(data$qc)) out <- c(out, sprintf("Genetic samples: %d of %d samples passed laboratory quality control and were used; the analysis assumes every retained genotype is assigned to the correct individual (no allelic dropout or false alleles creating 'ghost' individuals, no shadow effect). Identity errors are not modelled (see README).",
                                               data$qc$n_pass, data$qc$n_samples))
  out <- c(out,
    sprintf("Assumptions: the population is closed during sampling (no births, deaths or permanent movements), activity centres are fixed, detections of different individuals are independent, and detection declines with distance as a half-normal curve%s. Check the posterior predictive tests in Diagnostics.",
            if (prep$S == 2) " that differs between the sexes as modelled" else " shared by all individuals; strong sex differences in movement (common in carnivores) violate this, so fit sex-specific detection when sex is known"),
    paste0("Validation: at the module's default design (7x7 traps at 1.5 sigma spacing, 6 occasions, ~42 of ~87 animals detected) and in the Poisson hair-snare (effort matrix, habitat mask with a lake, two density covariates, trap-specific response), sex-specific (20% or 60% of animals of unknown sex) and Bernoulli (occasion covariate, partial effort) designs of VALIDATION.md, the estimates recovered the simulated truth within the acceptance tolerance (bias within 10%, ~95% interval coverage; Nimble, and NumPyro and secr where tested). ",
           "Posterior-mean density is biased upward by 12-20% when traps are 3 sigma apart with few spatial recaptures, or whenever density is imprecise (posterior CV of density above ~0.2: ~18-33 detected animals in the camera designs, ~46 in a sex-specific hair-snare design with two density covariates)",
           if (prep$obs != "poisson" && prep$behav == "trap") "; with 0/1 data and a trap-specific response, baseline detection p0 was biased upward by ~16% and the response toward 'trap-shy' (density, abundance and sigma were unaffected)" else "",
           ". See modules/scr/VALIDATION.md for these and the other documented limits."))
  out
}

# ---------------------------------------------------------------- figures
# Human label for any reported parameter
scr_param_label <- function(p, data, opts, prep) {
  ul <- scr_unit_len(prep); al <- scr_area_label(prep$units)
  base <- c(density_per_km2 = "Density (per km²)", density_per_unit2 = "Density (per map unit²)", abundance = "Abundance (N)",
            baseline_detection = "Baseline detection p0", baseline_encounter_rate = "Baseline encounter rate lambda0",
            baseline_detection_female = "p0, females", baseline_detection_male = "p0, males",
            baseline_encounter_rate_female = "lambda0, females", baseline_encounter_rate_male = "lambda0, males",
            sigma = sprintf("sigma (%s)", ul), sigma_female = sprintf("sigma, females (%s)", ul), sigma_male = sprintf("sigma, males (%s)", ul),
            home_range_area_95 = sprintf("95%% home range (%s)", al), home_range_area_95_female = sprintf("95%% home range, females (%s)", al),
            home_range_area_95_male = sprintf("95%% home range, males (%s)", al), proportion_female = "Proportion female",
            detection_beta_male = "Detection: male vs female", sigma_beta_male = "Sigma: male vs female",
            detection_beta_trap_response = "Detection: after first detection at the detector")
  vapply(p, function(x) {
    lab <- if (x %in% names(base)) base[[x]] else
      if (grepl("^density_beta_", x)) sprintf("Density: %s", scr_lab(sub("^density_beta_", "", x), data, opts)) else
      if (grepl("^sigma_beta_", x)) sprintf("Sigma: %s", scr_lab(sub("^sigma_beta_", "", x), data, opts)) else
      if (grepl("^detection_beta_", x)) sprintf("Detection: %s", scr_lab(sub("^detection_beta_", "", x), data, opts)) else x
    fig_label(x, lab, opts)
  }, "")
}

scr_band <- function(P, x, a, ai) data.frame(x = x, estimate = colMeans(P), lower = apply(P, 2, stats::quantile, a), upper = apply(P, 2, stats::quantile, 1 - a),
                                             lower50 = apply(P, 2, stats::quantile, ai), upper50 = apply(P, 2, stats::quantile, 1 - ai))
scr_iv <- function(v, a) c(estimate = mean(v), lower = stats::quantile(v, a, names = FALSE), upper = stats::quantile(v, 1 - a, names = FALSE))

scr_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  ob <- fit$object; prep <- ob$prep
  pois <- prep$obs == "poisson"
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci; a <- (1 - ci) / 2
  ici <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci; ai <- (1 - ici) / 2
  al <- scr_area_label(prep$units); ul <- scr_unit_len(prep)
  bandtype <- if (ob$engine == "secr") "confidence (parametric bootstrap)" else "credible"
  pd <- scr_pdraws(fit, 1000)
  figs <- list()
  # map coordinates back to user units
  # metre coordinates spanning > 5 km (e.g. UTM) are drawn in km so axis labels do not collide
  span <- max(diff(range(prep$mask$x)), diff(range(prep$mask$y)))
  mk <- if (identical(prep$units, "m") && span > 5000) 1000 else 1
  ulm <- if (mk == 1000) "km" else ul
  to_x <- function(v) (v * prep$cscale + prep$origin[1]) / mk; to_y <- function(v) (v * prep$cscale + prep$origin[2]) / mk
  xl <- fig_label("easting", sprintf("Easting (%s)", ulm), opts); yl <- fig_label("northing", sprintf("Northing (%s)", ulm), opts)
  map_theme <- function(p) p + ggplot2::scale_x_continuous(n.breaks = 4) + ggplot2::scale_y_continuous(n.breaks = 4) +
    ggplot2::theme(legend.position = "right")
  trap_pts <- data.frame(x = prep$traps$x / mk, y = prep$traps$y / mk)
  dword <- if (pois) "detectors" else "traps"
  surf <- ob$surfaces
  if (is.null(surf)) surf <- tryCatch(scr_surfaces(prep, pd[scr_thin_idx(nrow(pd), 100), , drop = FALSE] |> scr_with_M(pd)), error = function(e) NULL)
  if (!is.null(surf)) {
    cellA <- prep$pix^2 / scr_area_conv(prep)
    grid <- data.frame(x = prep$mask$x / mk, y = prep$mask$y / mk, value = surf$realised / cellA)
    p <- map_theme(fig_map(grid, sprintf("Animals per %s", al), points = trap_pts) + ggplot2::labs(x = xl, y = yl))
    if (isTRUE(opts$show_data) && prep$n) {
      # most probable habitat pixel of each activity centre (a posterior MEAN can fall in
      # non-habitat, e.g. a lake, when the posterior straddles it)
      cen <- prep$mxy[apply(surf$Pc, 1, which.max), , drop = FALSE]
      p <- p + ggplot2::geom_point(data = data.frame(x = to_x(cen[, 1]), y = to_y(cen[, 2])), ggplot2::aes(x, y),
                                   inherit.aes = FALSE, shape = 4, colour = "#D55E00", size = 1.1, stroke = 0.5)
    }
    figs$density_surface <- with_caption(p, sprintf(
      "Posterior mean realised density (animals per %s) over the habitat mask (non-habitat left blank): detected animals placed by their estimated activity centres plus the expected density of undetected animals. White circles: %s; orange crosses: most probable activity-centre pixel of each of the %d detected individuals.", al, dword, prep$n))
    attr(figs$density_surface, "fig_size") <- list(width = "double", height = 130)
    if (length(prep$dcov)) {
      g2 <- data.frame(x = prep$mask$x / mk, y = prep$mask$y / mk, value = surf$expected / cellA)
      figs$density_expected_surface <- with_caption(
        map_theme(fig_map(g2, sprintf("Expected animals per %s", al), points = trap_pts) + ggplot2::labs(x = xl, y = yl)),
        sprintf("Expected density surface implied by the habitat model (%s); white circles: %s.", paste(vapply(prep$dcov, scr_lab, "", data = data, opts = opts), collapse = ", "), dword))
      attr(figs$density_expected_surface, "fig_size") <- list(width = "double", height = 130)
    }
  }
  # conditional effects of density covariates (response scale)
  for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]; sc <- prep$dscale[[v]]
    xs <- seq(min(prep$mask[[v]]), max(prep$mask[[v]]), length.out = 60)
    Dm <- sapply(xs, function(x) { z <- rep(0, length(prep$dcov)); z[k] <- (x - sc$center) / sc$scale; scr_density_at(prep, pd, z) })
    obs <- if (isTRUE(opts$show_data)) data.frame(x = prep$mask[[v]][scr_thin_idx(prep$G, 300)]) else NULL
    figs[[paste0("density_", v)]] <- with_caption(
      fig_cond_effect(scr_band(Dm, xs, a, ai), scr_lab(v, data, opts), fig_label("density", sprintf("Density (animals per %s)", al), opts), obs = obs) +
        ggplot2::scale_y_continuous(limits = c(0, NA)),
      sprintf("Expected density versus %s (other density covariates at their habitat mean), with %d%% and %d%% %s bands. Rug: covariate values across the habitat (the range the model can speak to).",
              scr_lab(v, data, opts), round(100 * ci), round(100 * ici), bandtype))
  }
  # ---- detection function (by sex, or by the first 0/1 detector covariate)
  b <- prep$bidx
  link <- if (pois) log else stats::qlogis; ilink <- if (pois) exp else stats::plogis
  bcol <- function(i) if (i > 0) pd[, paste0("beta", i)] else 0
  dmax <- 4 * stats::median(pd[, "sigma"] * if (prep$S == 2) exp(bcol(b[8])) else 1) * prep$cscale
  dd <- seq(0, dmax, length.out = 80)
  bin <- prep$pcov[vapply(prep$pcov, function(v) isTRUE(prep$pscale[[v]]$binary), TRUE)]
  groups <- if (prep$S == 2) list(Female = list(e = 0, s = 0), Male = list(e = bcol(b[7]), s = bcol(b[8]))) else
    if (length(bin)) { bi <- b[3] + match(bin[1], c(prep$pcov, prep$ocov))
      stats::setNames(list(list(e = 0, s = 0), list(e = pd[, paste0("beta", bi)], s = 0)), sprintf("%s = %d", scr_lab(bin[1], data, opts), 0:1)) } else
    list(all = list(e = 0, s = 0))
  rows <- list()
  for (gn in names(groups)) {
    lp <- link(pd[, "det0"]) + groups[[gn]]$e; sg <- pd[, "sigma"] * exp(groups[[gn]]$s)
    P <- sapply(dd, function(x) ilink(lp) * exp(-(x / prep$cscale)^2 / (2 * sg^2)))
    r <- scr_band(P, dd, a, ai)
    if (length(groups) > 1) r$group <- fig_label(gn, gn, opts)
    rows[[length(rows) + 1]] <- r
  }
  ylab_det <- if (pois) fig_label("encounter_rate", "Expected samples per occasion", opts) else fig_label("detection_prob", "Detection prob. per occasion", opts)
  pdet <- fig_cond_effect(do.call(rbind, rows), fig_label("distance", sprintf("Distance from activity centre to detector (%s)", ul), opts), ylab_det) +
    ggplot2::geom_vline(xintercept = prep$trap_spacing, linetype = 2, colour = "grey45") +
    ggplot2::annotate("text", x = prep$trap_spacing, y = Inf, label = " detector spacing", hjust = 0, vjust = 1.3, size = 2.6, colour = "grey30") +
    ggplot2::scale_y_continuous(limits = c(0, NA))
  if (length(groups) > 1) pdet <- pdet + ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL) + ggplot2::theme(legend.position = "top")
  attr(pdet, "fig_size") <- list(aspect = 0.8)
  figs$detection_function <- with_caption(pdet, sprintf(
    "Half-normal detection function: %s at distance d, %s exp(-d^2 / 2 sigma^2), with %d%% and %d%% %s bands%s%s. Dashed line: detector spacing; spacing should be <= 2 sigma (Sun et al. 2014).",
    if (pois) "expected number of samples per occasion (per unit effort)" else "probability of detection on one occasion",
    if (pois) "lambda0" else "p0", round(100 * ci), round(100 * ici), bandtype,
    if (length(groups) > 1) sprintf(", by %s", if (prep$S == 2) "sex" else scr_lab(bin[1], data, opts)) else "",
    if (length(c(prep$pcov, prep$ocov, prep$scov))) "; other covariates at their reference level" else ""))
  # ---- detection covariates (baseline detection at d = 0, response scale)
  pc <- c(prep$pcov, prep$ocov)
  for (q in seq_along(pc)) {
    v <- pc[q]; isocc <- v %in% prep$ocov
    sc <- if (isocc) prep$oscale[[v]] else prep$pscale[[v]]
    raw <- if (isocc) as.numeric(data$occasion_covariates[[v]]) else as.numeric(prep$traps[[v]])
    bq <- pd[, paste0("beta", b[3] + q)]
    ylab <- if (pois) fig_label("baseline_encounter_rate", "Baseline encounter rate lambda0", opts) else fig_label("baseline_detection", "Baseline detection p0", opts)
    if (isTRUE(sc$binary)) {
      pr <- do.call(rbind, lapply(0:1, function(l) data.frame(level = sprintf("%d", l), t(scr_iv(ilink(link(pd[, "det0"]) + l * bq), a)))))
      pr$level <- factor(pr$level, levels = c("0", "1"), labels = c(fig_label(paste0(v, "_0"), "No (0)", opts), fig_label(paste0(v, "_1"), "Yes (1)", opts)))
      pp <- fig_cat_effect(pr, scr_lab(v, data, opts), ylab) + ggplot2::scale_y_continuous(limits = c(0, NA))
    } else {
      xs <- seq(min(raw), max(raw), length.out = 50)
      P <- sapply(xs, function(x) ilink(link(pd[, "det0"]) + (x - sc$center) / sc$scale * bq))
      pp <- fig_cond_effect(scr_band(P, xs, a, ai), scr_lab(v, data, opts), ylab, obs = if (isTRUE(opts$show_data)) data.frame(x = raw) else NULL) +
        ggplot2::scale_y_continuous(limits = c(0, NA))
    }
    figs[[paste0("detection_", v)]] <- with_caption(pp, sprintf(
      "%s versus %s (%s-level covariate; %s), with %d%% %s %s%s. Detection effects describe how easily animals are recorded, not where they live.",
      if (pois) "Baseline encounter rate (expected samples per occasion at distance 0)" else "Baseline detection probability per occasion at distance 0",
      scr_lab(v, data, opts), if (isocc) "occasion" else "detector", if (prep$S == 2) "females; other covariates at reference" else "other covariates at their reference level",
      round(100 * ci), bandtype, if (isTRUE(sc$binary)) "intervals" else "bands",
      if (!isTRUE(sc$binary) && isTRUE(opts$show_data)) sprintf("; rug: observed %s values", if (isocc) "occasion" else "detector") else ""))
  }
  # ---- behavioural response
  if (b[9] > 0) {
    bb <- pd[, paste0("beta", b[9])]
    pr <- rbind(data.frame(level = "Not yet detected there", t(scr_iv(ilink(link(pd[, "det0"])), a))),
                data.frame(level = "After first detection there", t(scr_iv(ilink(link(pd[, "det0"]) + bb), a))))
    pr$level <- factor(pr$level, levels = pr$level)
    figs$detection_trap_response <- with_caption(
      fig_cat_effect(pr, NULL, if (pois) "Baseline encounter rate lambda0" else "Baseline detection p0") + ggplot2::scale_y_continuous(limits = c(0, NA)),
      sprintf("Trap-specific behavioural response: baseline %s at a detector before and after the animal's first detection at that detector (%d%% %s intervals). Higher after = 'trap-happy' (e.g. bait, scent), lower = 'trap-shy'.",
              if (pois) "encounter rate" else "detection probability", round(100 * ci), bandtype))
  }
  # ---- sex differences
  if (prep$S == 2) {
    mk_sex <- function(f, m, ylab) {
      pr <- rbind(data.frame(level = "Female", t(scr_iv(f, a))), data.frame(level = "Male", t(scr_iv(m, a))))
      pr$level <- factor(pr$level, levels = c("Female", "Male"), labels = c(fig_label("Female", "Female", opts), fig_label("Male", "Male", opts)))
      fig_cat_effect(pr, NULL, ylab) + ggplot2::scale_y_continuous(limits = c(0, NA))
    }
    d0 <- pd[, "det0"]
    pl <- list(mk_sex(d0, ilink(link(d0) + bcol(b[7])), if (pois) "lambda0" else "p0"),
               mk_sex(pd[, "sigma"] * prep$cscale, pd[, "sigma"] * exp(bcol(b[8])) * prep$cscale, sprintf("sigma (%s)", ul)),
               ggplot2::ggplot(data.frame(phi = pd[, "phi"]), ggplot2::aes(phi)) +
                 ggplot2::geom_histogram(bins = 30, fill = eco_pal$band, colour = "white", alpha = 0.7) +
                 ggplot2::geom_vline(xintercept = 0.5, linetype = 2, colour = eco_pal$null) +
                 ggplot2::labs(x = "Proportion female", y = "Draws"))
    bs <- if (is.null(opts$base_size)) 9 else opts$base_size
    pl <- lapply(pl, function(p) p + theme_eco(base_size = bs))   # export themes only the last panel of a patchwork
    figs$sex_differences <- with_caption(compose_figures(pl, ncol = 3),
      sprintf("Sex differences: (A) baseline %s, (B) spatial scale sigma, by sex (points: %s; bars: %d%% %s intervals); (C) distribution of the proportion of females in the population (dashed: 1:1 sex ratio). Animals of unknown sex (%d of %d detected) contribute through a mixture over the two sexes.",
              if (pois) "encounter rate lambda0" else "detection p0", if (ob$engine == "secr") "MLE-based" else "posterior means", round(100 * ci), bandtype,
              sum(prep$sexi < 0), prep$n))
    attr(figs$sex_differences, "fig_size") <- list(width = "double", aspect = 0.35)
  }
  # ---- detector effort / samples (NGS)
  if (pois || !prep$binary_usage) {
    ef <- data.frame(x = prep$traps$x / mk, y = prep$traps$y / mk, effort = rowSums(prep$usage), n = colSums(prep$y))
    figs$detector_effort <- with_caption(
      ggplot2::ggplot(ef, ggplot2::aes(x, y)) +
        ggplot2::geom_point(ggplot2::aes(size = effort, fill = n), shape = 21, colour = "grey20", stroke = 0.3) +
        scale_fill_eco_seq() + ggplot2::scale_size(range = c(0.8, 4.5)) + ggplot2::coord_equal() +
        ggplot2::scale_x_continuous(n.breaks = 4) + ggplot2::scale_y_continuous(n.breaks = 4) +
        ggplot2::labs(x = xl, y = yl, size = "Total effort", fill = if (pois) "Samples" else "Detections"),
      sprintf("Detector effort and %s: symbol size = total effort summed over the %d occasions (1 = one full occasion; occasion effort is 0 when a detector was not operating), colour = number of %s. Effort enters the model as a multiplier of the detection hazard (Efford et al. 2013).",
              if (pois) "samples" else "detections", prep$K, if (pois) "genetic samples (after QC)" else "detections"))
    asp <- diff(range(ef$y)) / max(diff(range(ef$x)), 1e-9)
    attr(figs$detector_effort, "fig_size") <- list(width = "onehalf", aspect = min(1.1, max(0.5, 0.15 + 0.8 * asp)))
  }
  # ---- abundance
  if (!is.null(ob$rep_chains)) {
    Nd <- unlist(lapply(ob$rep_chains, function(m) m[, "abundance"]))
    qq <- scr_q(Nd, c(a, 1 - a))
    figs$abundance_posterior <- with_caption(
      ggplot2::ggplot(data.frame(N = Nd), ggplot2::aes(N)) +
        ggplot2::geom_histogram(binwidth = max(1, round(diff(range(Nd)) / 40)), fill = eco_pal$band, colour = "white", alpha = 0.7) +
        ggplot2::geom_vline(xintercept = prep$n, linetype = 3, colour = "grey30") +
        ggplot2::geom_vline(xintercept = qq, linetype = 2, colour = eco_pal$ink) +
        ggplot2::labs(x = fig_label("abundance", "Abundance in the state space (N)", opts), y = "Posterior draws"),
      sprintf("Posterior distribution of the number of activity centres in the %s %s of habitat in the state space; dashed lines: %d%% credible interval; dotted line: number of individuals detected (%d).",
              scr_fmt(prep$area * scr_hr_conv(prep), 3), al, round(100 * ci), prep$n))
  } else {
    r <- fit$estimates[fit$estimates$param == "abundance", ]
    figs$abundance_posterior <- with_caption(
      ggplot2::ggplot(r, ggplot2::aes(estimate, "N")) + ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper)) +
        ggplot2::geom_vline(xintercept = prep$n, linetype = 3) + ggplot2::labs(x = "Abundance in the state space (N)", y = NULL),
      "Realised abundance in the state space (secr region.N) with 95% confidence interval; dotted line: individuals detected.")
  }
  # ---- captures per individual (observed vs predictive)
  ppc <- ob$ppc
  if (!is.null(ppc) && prep$n) {
    obsn <- table(factor(pmin(rowSums(prep$y), 10), levels = 1:10))
    repn <- t(vapply(ppc$caps_rep, function(v) as.numeric(table(factor(pmin(v, 10), levels = 1:10))), numeric(10)))
    dfc <- data.frame(k = factor(c(1:9, "10+"), levels = c(1:9, "10+")), obs = as.numeric(obsn),
                      est = colMeans(repn), lo = apply(repn, 2, stats::quantile, a), hi = apply(repn, 2, stats::quantile, 1 - a))
    figs$captures_per_individual <- with_caption(
      ggplot2::ggplot(dfc, ggplot2::aes(k)) +
        ggplot2::geom_col(ggplot2::aes(y = obs), fill = "grey75", width = 0.7) +
        ggplot2::geom_pointrange(ggplot2::aes(y = est, ymin = lo, ymax = hi), colour = eco_pal$ink, size = 0.3) +
        ggplot2::labs(x = if (pois) "Samples per detected individual" else "Detections per detected individual", y = "Number of individuals"),
      sprintf("Bars: observed number of individuals %s 1, 2, ... times. Points and lines: %s predictive mean and %d%% interval from data simulated under the fitted model. Bars far outside the intervals indicate unmodelled individual heterogeneity in detection.",
              if (pois) "sampled" else "detected", if (ob$engine == "secr") "parametric" else "posterior", round(100 * ci)))
    nm <- c(n_detected = "Individuals detected", total_detections = if (pois) "Total samples" else "Total detections",
            traps_per_individual = "Detectors per individual", disp_detections = "Dispersion (var/mean) per individual")
    pl <- lapply(names(nm), function(s) fig_ppc(ppc$rep[, s], stats::median(ppc$obs[, s], na.rm = TRUE), nm[[s]]))
    ftd <- data.frame(o = ppc$obs[, "ft_trap_totals"], r = ppc$rep[, "ft_trap_totals"])
    pl[[5]] <- ggplot2::ggplot(ftd, ggplot2::aes(o, r)) + ggplot2::geom_point(alpha = 0.4, size = 0.8, colour = eco_pal$ink) +
      ggplot2::geom_abline(linetype = 2, colour = "grey45") +
      ggplot2::labs(x = "Freeman-Tukey of detector totals, observed", y = "FT, replicated",
                    caption = sprintf("P(replicate >= observed) = %.2f", ppc$p[["ft_trap_totals"]]))
    pl <- lapply(pl, function(p) p + theme_eco(base_size = 8))
    figs$gof <- with_caption(compose_figures(pl, ncol = 2),
      "Predictive checks: histograms of statistics computed from data simulated under the fitted model, line = observed value; last panel: Freeman-Tukey discrepancy of detector totals, observed vs replicated (points above the 1:1 line = data consistent with the model).")
    attr(figs$gof, "fig_size") <- list(width = "double", aspect = 1)
  }
  # ---- coefficients (link scale)
  ec <- fit$estimates[grepl("_beta_", fit$estimates$param), ]
  if (nrow(ec)) {
    ec$param <- scr_param_label(ec$param, data, opts, prep)
    pco <- fig_coef(ec, null_value = 0, xlab = sprintf("Coefficient (log for density and sigma, %s for detection)", if (pois) "log" else "logit"), opts = opts) +
      ggplot2::scale_y_discrete(labels = function(x) vapply(strwrap(x, 28, simplify = FALSE), paste, "", collapse = "\n"))
    attr(pco, "fig_size") <- list(height = 45 + 12 * nrow(ec))
    figs$coefficients <- with_caption(pco, sprintf(
      "Link-scale coefficients with 95%% intervals: density covariates on the log scale per 1 SD of the covariate; detection covariates on the %s scale (0/1 covariates: level 1 vs 0; continuous: per 1 SD); sigma covariates on the log scale; sex effects are male vs female. Dashed line: no effect.", if (pois) "log" else "logit"))
  }
  # ---- buffer check
  if (!is.null(ob$buffer_check)) {
    bc <- ob$buffer_check
    bc$lab <- sprintf("%s\n(%.1f sigma)", bc$fit, bc$buffer_sigma)
    bc$lab <- factor(bc$lab, levels = bc$lab[order(bc$buffer)])
    figs$buffer_check <- with_caption(
      ggplot2::ggplot(bc, ggplot2::aes(lab, estimate)) +
        ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), colour = eco_pal$ink, size = 0.35) +
        ggplot2::scale_y_continuous(limits = c(0, NA)) +
        ggplot2::labs(x = "State-space buffer beyond the detectors", y = fig_label("density", sprintf("Density (animals per %s)", al), opts)),
      sprintf("Buffer check: density in the region within %s %s of the detectors, estimated with the main state space and with a state space whose buffer was changed by 2 sigma (95%% intervals). Similar estimates show the buffer is wide enough; density that rises as the buffer shrinks shows it is too small.",
              if (is.null(bc$region_buffer)) "the smaller buffer" else format(signif(bc$region_buffer[1], 3)), if (is.null(bc$region_buffer)) "" else ul))
  }
  # ---- engine comparison
  if (!is.null(ob$comparison)) {
    cm <- ob$comparison
    cm$param <- scr_param_label(cm$param, data, opts, prep)
    eng_desc <- c("nimble marginal" = "Bayesian marginal likelihood (Nimble)", "nimble latent" = "Bayesian latent data augmentation (Nimble)",
                  "numpyro" = "Bayesian marginal NUTS (NumPyro)", "secr" = "maximum likelihood (secr)")
    ue <- unique(cm$engine); ed <- ifelse(ue %in% names(eng_desc), eng_desc[ue], ue)
    figs$engine_comparison <- with_caption(
      ggplot2::ggplot(cm, ggplot2::aes(estimate, engine, colour = engine, shape = engine)) +
        ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.25) +
        ggplot2::facet_wrap(~param, scales = "free_x", labeller = ggplot2::label_wrap_gen(22)) +
        ggplot2::scale_x_continuous(breaks = function(l) { b <- pretty(l, 3); b[b >= l[1] + 0.1 * diff(l) & b <= l[2] - 0.1 * diff(l)] }) +
        ggplot2::labs(x = "Estimate (95% interval)", y = NULL) +
        ggplot2::theme(panel.spacing.x = ggplot2::unit(1.5, "lines")) +
        ggplot2::guides(colour = "none", shape = "none") + ggplot2::theme(legend.position = "none"),
      sprintf("The same data fitted with %d engines: %s. Points: posterior mean (Bayesian) or MLE (secr); bars: 95%% interval. Agreement indicates the result does not depend on the software.",
              length(ue), paste(ed, collapse = "; ")))
    attr(figs$engine_comparison, "fig_size") <- list(width = "double", aspect = 0.8)
  }
  if (!is.null(ob$rep_chains)) {
    figs$trace <- with_caption(fig_trace(ob$rep_chains), "MCMC traces of the reported quantities, one colour per chain (thinned).")
    attr(figs$trace, "fig_size") <- list(width = "double", aspect = 0.8)
    figs$posterior <- with_caption(fig_posterior(do.call(rbind, ob$rep_chains)), "Posterior densities of the reported quantities.")
    attr(figs$posterior, "fig_size") <- list(width = "double", aspect = 0.8)
  }
  figs
}

scr_with_M <- function(x, pd) { attr(x, "M") <- attr(pd, "M"); x }

# ---------------------------------------------------------------- methods
scr_methods_text <- function(fit, data, settings) {
  ob <- fit$object; prep <- ob$prep; st <- ob$settings
  al <- scr_area_label(prep$units); pois <- prep$obs == "poisson"
  dens <- if (length(prep$dcov)) sprintf("log-linear in %s (standardised to mean 0, SD 1 over the habitat mask)", paste(prep$dcov, collapse = ", ")) else "homogeneous"
  obs <- switch(prep$obs,
    binomial = "The number of occasions on which individual i was detected at detector j was Binomial(K_j, p_ij), where K_j is the number of occasions detector j operated and p_ij",
    bernoulli = "Detection of individual i at detector j on occasion k was Bernoulli with probability 1 - (1 - p_ijk)^u_jk, where u_jk is the effort of detector j on occasion k (0 when not operating; effort acts on the hazard scale, Efford et al. 2013) and p_ijk",
    poisson = "The number of samples of individual i at detector j on occasion k was Poisson with mean u_jk lambda_ijk, where u_jk is the effort of detector j on occasion k (0 when not operating; Efford et al. 2013) and lambda_ijk")
  fn <- if (pois) "lambda0 exp(-d_ij^2 / (2 sigma^2)) is a hazard half-normal function (expected samples per occasion)" else "p0 exp(-d_ij^2 / (2 sigma^2)) is a half-normal function"
  extra <- character()
  pc <- c(prep$pcov, prep$ocov)
  if (length(pc)) extra <- c(extra, sprintf("%s(%s) varied with %s", if (pois) "log" else "logit", if (pois) "lambda0" else "p0",
                                            paste(c(if (length(prep$pcov)) sprintf("detector-level covariate(s) %s", paste(prep$pcov, collapse = ", ")),
                                                    if (length(prep$ocov)) sprintf("occasion-level covariate(s) %s", paste(prep$ocov, collapse = ", "))), collapse = " and ")))
  if (length(prep$scov)) extra <- c(extra, sprintf("log(sigma) varied with detector-level covariate(s) %s", paste(prep$scov, collapse = ", ")))
  if (prep$behav == "trap") extra <- c(extra, sprintf("a trap-specific behavioural response shifted %s(%s) at a detector after an individual's first detection there", if (pois) "log" else "logit", if (pois) "lambda0" else "p0"))
  if (prep$S == 2) extra <- c(extra, sprintf("%s differed between the sexes; the proportion of females was estimated and individuals of unknown sex (%d of %d) were included by summing the likelihood over both sexes (sex as a partially observed individual covariate; Royle et al. 2014)",
                                             paste(c(if ("p0" %in% prep$sexeff) if (pois) "lambda0" else "p0", if ("sigma" %in% prep$sexeff) "sigma"), collapse = " and "), sum(prep$sexi < 0), prep$n))
  base <- sprintf(paste(
    "We estimated population density with a spatial capture-recapture (SCR) model (Efford 2004; Borchers & Efford 2008; Royle et al. 2014).",
    "Each animal has a fixed activity centre in a state space represented by a habitat mask of %d pixels (%s %s spacing, %s %s of habitat, extending %s %s beyond the detectors%s).",
    "Activity centres follow a Poisson point process whose intensity is %s.",
    "%s = %s of the distance from the activity centre to the detector%s.",
    "The 95%% home-range area was derived as pi (2.448 sigma)^2 assuming circular bivariate-normal space use."),
    prep$G, format(signif(prep$spacing, 3)), scr_unit_len(prep), scr_fmt(prep$area * scr_hr_conv(prep), 3), al,
    if (is.null(prep$buffer)) "?" else format(signif(prep$buffer, 3)), scr_unit_len(prep),
    switch(if (is.null(data$mask_kind)) "csv" else data$mask_kind, polygon = "; pixels outside the habitat polygon excluded", raster = "; non-habitat raster cells excluded", ""),
    dens, obs, fn, if (length(extra)) paste0("; ", paste(extra, collapse = "; ")) else "")
  qc <- if (!is.null(data$qc)) sprintf(" Only genetic samples that passed laboratory quality control were analysed (%d of %d).", data$qc$n_pass, data$qc$n_samples) else ""
  eng <- if (ob$engine == "secr") sprintf("The model was fitted by maximum likelihood in the R package secr (version %s; Efford 2025) with %s detectors and a %s detection function%s; realised abundance in the state space and its interval were obtained with region.N.",
                                          as.character(utils::packageVersion("secr")), if (pois) "count (Poisson)" else "proximity", if (pois) "hazard half-normal" else "half-normal",
                                          if (prep$S == 2) ", with sex as a hybrid-mixture class (hcov)" else "")
  else sprintf(paste(
    "The model was fitted in a Bayesian framework with data augmentation (Royle & Dorazio 2012; M = %d potential individuals, inclusion probability psi ~ Uniform(0,1)) with priors %s, sigma ~ Uniform(0, %s %s)%s and Normal(0, SD 2) for all coefficients.",
    "%s",
    "Abundance N was derived per draw as n + Binomial(M - n, P(z = 1 | never detected)) and density as psi M / area.",
    "Convergence was assessed with rank-normalised R-hat and bulk effective sample size (posterior package); model fit was checked with posterior predictive tests of the number of individuals detected, total detections, detectors per individual, dispersion of detections per individual and a Freeman-Tukey statistic of detector totals."),
    prep$M_used, if (pois) sprintf("lambda0 ~ Uniform(0, %s)", format(prep$lam_max)) else "p0 ~ Uniform(0,1)",
    format(signif(prep$sigma_max * prep$cscale, 3)), scr_unit_len(prep), if (prep$S == 2) ", proportion female ~ Uniform(0,1)" else "",
    switch(sub(" .*", "", ob$engine),
      nimble = if (grepl("latent", ob$engine)) sprintf("Activity centres and inclusion indicators were sampled explicitly (latent data augmentation) in NIMBLE %s (de Valpine et al. 2017), %d chains of %d iterations after %d burn-in, thinned by %d.",
                                                       as.character(utils::packageVersion("nimble")), st$nchains, st$niter_latent - st$nburnin_latent, st$nburnin_latent, st$thin_latent)
               else sprintf("Activity centres%s and inclusion indicators were integrated out over the habitat mask (marginal likelihood) and the remaining parameters sampled jointly with an %s sampler in NIMBLE %s (de Valpine et al. 2017), %d chains of %d iterations after %d burn-in.",
                            if (prep$S == 2) ", unknown sexes" else "",
                            if (identical(st$sampler, "AF_slice")) "automated-factor slice" else "adaptive block random-walk", as.character(utils::packageVersion("nimble")), st$nchains, st$niter - st$nburnin, st$nburnin),
      numpyro = sprintf("Activity centres and inclusion indicators were integrated out over the habitat mask and the parameters sampled with the No-U-Turn sampler in NumPyro (Phan et al. 2019), %d chains of %d draws after %d warm-up.",
                        st$nchains, st$numpyro_samples, st$numpyro_warmup)))
  bc <- if (!is.null(ob$buffer_check)) sprintf(" Adequacy of the state-space buffer was checked by refitting with a buffer of %.1f instead of %.1f sigma (density changed by %+.1f%%).",
                                               ob$buffer_check$buffer_sigma[2], ob$buffer_check$buffer_sigma[1], 100 * (ob$buffer_check$estimate[2] / ob$buffer_check$estimate[1] - 1)) else ""
  paste0(base, qc, " ", eng, bc)
}
