# jolly-seber: open-population capture-recapture (Jolly 1965; Seber 1965;
# Crosbie & Manly 1985; Schwarz & Arnason 1996 POPAN superpopulation model)
# fitted by data augmentation in the restricted dynamic occupancy
# formulation (Royle & Dorazio 2008, 2012; Kery & Schaub 2012 ch. 10).
#
# Abundance per occasion N_t, superpopulation size, entry probabilities,
# entrants and per-capita recruitment, apparent survival (optionally by group
# such as sex, or by interval), capture probability (constant, by occasion or
# with occasion covariates such as effort) and realised growth rates.
# Engines (settings$engine):
#   nimble  approach "marginal" (default; latent states summed out by the
#           forward algorithm, N_t etc. drawn exactly afterwards) or
#           "latent" (every z[i,t] sampled; the textbook RDO code)
#   marked  maximum-likelihood POPAN (marked::crm model "JS"), cross-check

for (f in c("simulate.R", "prepare.R", "nimble_engines.R", "post.R", "marked_engine.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

js_settings <- list(
  engine = "nimble", approach = "marginal", survival = "constant", detection = "auto",
  group_effect = TRUE, M = NULL, prior_sd_slope = 2, seed = 1,
  niter = 20000, nburnin = 5000, nchains = 3, thin = 5,
  niter_latent = 20000, nburnin_latent = 5000, thin_latent = 5,
  max_draws_per_chain = 1500, ppc_draws = 300, gof = TRUE, compare_engines = character())

js_fast <- merge_list(js_settings, list(niter = 6000, nburnin = 2000, thin = 2, max_draws_per_chain = 500,
                                        ppc_draws = 0, gof = FALSE, niter_latent = 8000, nburnin_latent = 3000,
                                        thin_latent = 2))

js_bayes <- function(prep, settings) {
  for (attempt in 1:3) {
    eng <- js_engine_nimble(prep, settings)
    mx <- settings$max_draws_per_chain
    rep_chains <- lapply(eng$chains, function(m) {
      S <- m[js_thin_idx(nrow(m), mx), , drop = FALSE]
      lat <- NULL
      if (eng$approach == "latent") {
        T <- prep$T
        lat <- list(N = S[, paste0("N[", seq_len(T), "]"), drop = FALSE],
                    B = S[, paste0("B[", seq_len(T - 1), "]"), drop = FALSE], Nsuper = S[, "Nsuper"],
                    Nsuper2 = S[, "Nsuper2"])
      }
      js_report_chain(prep, S, lat)
    })
    Nq <- stats::quantile(unlist(lapply(rep_chains, function(m) m[, "superpopulation_size"])), 0.995)
    if (Nq < 0.9 * prep$M || attempt == 3) break
    newM <- 2L * prep$M
    prep$notes <- c(prep$notes, sprintf("Data augmentation limit M = %d was too small (99.5%% quantile of the superpopulation = %.0f); refitted with M = %d.",
                                        prep$M, Nq, newM))
    prep <- js_set_M(prep, newM)
  }
  s <- summarise_draws_eco(rep_chains)
  list(eng = eng, prep = prep, rep_chains = rep_chains, estimates = s$estimates, diag = s$diagnostics,
       M_hit = Nq >= 0.9 * prep$M, Nq = Nq)
}

js_fit <- function(data, settings) {
  settings <- merge_list(js_settings, settings)
  set.seed(settings$seed)
  prep <- js_prepare(data, settings)
  if (prep$n < 2) stop("fewer than 2 animals were captured: Jolly-Seber cannot be fitted")
  if (settings$engine == "marked") {
    e <- js_engine_marked(prep, settings)
    est <- e$estimates
    dg <- list(converged = e$converged, notes = c(prep$notes, e$notes), t_total = e$t_total)
    obj <- list(engine = "marked (ML POPAN)", prep = prep, marked = e$fit, settings = settings)
  } else {
    b <- js_bayes(prep, settings)
    prep <- b$prep
    est <- b$estimates
    dg <- b$diag; dg$t_total <- b$eng$t_total
    dg$notes <- prep$notes
    dg$notes <- c(dg$notes, sprintf("Data augmentation: M = %d potential animals; 99.5%% posterior quantile of the superpopulation = %.0f (%s).",
                                    prep$M, b$Nq, if (b$M_hit) "TOO CLOSE to M" else "M is not limiting"))
    if (b$M_hit) { dg$converged <- FALSE
      dg$notes <- c(dg$notes, "The superpopulation still presses against the data-augmentation limit M: set settings$M larger.") }
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                      dg$rhat_max, dg$ess_min))
    obj <- list(engine = paste("nimble", b$eng$approach), prep = prep, chains = b$eng$chains,
                rep_chains = b$rep_chains, settings = settings)
    if (settings$ppc_draws > 0) {
      pool <- js_pool(b$eng$chains, settings$ppc_draws)
      obj$ppc <- js_ppc(prep, pool)
      dg$notes <- c(dg$notes, sprintf(
        "Posterior predictive checks (%d draws; Bayesian p-value = P(replicate discrepancy >= observed), extreme if < 0.025 or > 0.975): captures per animal %.2f, animals caught per occasion %.2f.",
        nrow(pool), obj$ppc$p[["individual"]], obj$ppc$p[["occasion"]]))
    }
  }
  if (isTRUE(settings$gof)) {
    obj$gof <- js_ucare(prep)
    if (!is.null(obj$gof)) dg$notes <- c(dg$notes, js_ucare_text(obj$gof))
  }
  dg$notes <- c(dg$notes, js_design_notes(prep, est))
  if (length(settings$compare_engines)) {
    cmp <- list(data.frame(engine = obj$engine, est[, c("param", "estimate", "sd", "lower", "upper")]))
    for (e2 in setdiff(unlist(settings$compare_engines), obj$engine)) {
      s2 <- settings; s2$compare_engines <- character(); s2$ppc_draws <- 0; s2$gof <- FALSE
      if (e2 %in% c("latent", "nimble latent")) { s2$engine <- "nimble"; s2$approach <- "latent"; s2$M <- prep$M }
      else if (e2 %in% c("marginal", "nimble marginal")) { s2$engine <- "nimble"; s2$approach <- "marginal" }
      else s2$engine <- e2
      r2 <- tryCatch(js_fit(data, s2), error = function(e) { message("compare engine ", e2, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(r2)) cmp[[length(cmp) + 1]] <- data.frame(engine = r2$object$engine, r2$estimates[, c("param", "estimate", "sd", "lower", "upper")])
    }
    obj$comparison <- do.call(rbind, cmp)
  }
  list(estimates = est, diagnostics = dg, object = obj)
}

js_defaults <- list(
  params = list(superpopulation_size = 350, entry_prob_first = 0.35, survival = 0.75,
                survival_beta_male = -0.5, prop_male = 0.5, detection = 0.45, detection_beta_effort = 0.5),
  design = list(n_occasions = 8, interval = 1, groups = TRUE, group_levels = c("female", "male"),
                group_name = "sex", prop_group_unknown = 0, detection_covariate = TRUE,
                effort_mean = 80, effort_sd = 20, temp_emigration = c(0, 0), first_year = 2017,
                time_unit = "year"))

module <- ecomodule(
  id = "jolly-seber",
  title = "Jolly-Seber open-population capture-recapture (abundance, survival, recruitment)",
  description = paste(
    "How many animals were present on each sampling occasion, how many used the site over the whole study,",
    "how many joined the population (births plus immigrants) between occasions, what share survived (and stayed),",
    "and is the population growing or declining? The Jolly-Seber model uses repeated captures of individually",
    "marked animals (tags, PIT tags, photo-ID, genetic ID) over several occasions (usually years or seasons) of an",
    "open population, and estimates abundance per occasion, superpopulation size, entry probabilities and per-capita",
    "recruitment, apparent survival (optionally by sex or group), capture probability (optionally with sampling effort)",
    "and realised population growth rates, accounting for animals that were never caught."),
  defaults = js_defaults,
  simulate = function(params, design, seed) { s <- js_simulate(params, design, seed); s[c("data", "truth")] },
  fit = js_fit, settings = js_settings, fast_settings = js_fast,
  read_data = js_read_data, check_data = js_check_data, write_data = js_write_data,
  interpret = js_interpret, figures = js_figures, plots = function(fit, data) js_figures(fit, data, fig_opts()),
  methods_text = js_methods_text,
  power_params = c("survival_beta_male", "detection_beta_effort"),
  engine = "nimble (data augmentation, RDO formulation; marginal [default] or latent), marked (ML POPAN cross-check)",
  references = js_references)
