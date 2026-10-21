# multistate: multistate (Arnason-Schwarz) and multievent (Pradel 2005)
# capture-recapture. State-specific survival, transitions between states
# (multinomial logit, "stay" as reference), state-specific detection, and,
# for multievent data, the probability that a detected animal's state is
# recorded plus the state distribution at first capture. The latent state
# sequence is summed out with the forward algorithm (nimbleEcology dDHMMo),
# conditional on first capture.

for (f in c("simulate.R", "prep.R", "engine_nimble.R", "loglik_r.R", "derive.R", "gof.R", "io.R", "output.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

ms_settings <- list(engine = "nimble", unknown_handling = "model", naive_state = NULL, model = NULL,
                    niter = 6000, nburnin = 2000, nchains = 3, thin = 1, prior_sd = 2,
                    sampler = "slice", seed = 1, ppc_draws = 300, ucare = TRUE, n_draws = 2000)
ms_fast <- merge_list(ms_settings, list(niter = 2000, nburnin = 700, ppc_draws = 0, ucare = FALSE))

ms_fit <- function(data, settings) {
  settings <- merge_list(ms_settings, settings)
  set.seed(settings$seed)
  des <- ms_prepare(data, settings)
  notes <- character()
  if (des$handling == "drop")
    notes <- c(notes, sprintf("%d 'seen, state unknown' records were treated as not seen (unknown_handling: drop).", des$n_unknown))
  if (des$handling == "assign")
    notes <- c(notes, sprintf("NEGATIVE CONTROL: %d 'seen, state unknown' records were forced into state '%s'. Transitions are expected to be biased.",
                              des$n_unknown, if (is.null(settings$naive_state)) des$states[1] else settings$naive_state))
  engine <- match.arg(settings$engine, c("nimble", "ml"))
  if (engine == "nimble") {
    eng <- ms_engine_nimble(des, settings)
    der <- lapply(eng$chains, ms_derive, des = des)
    s <- summarise_draws_eco(der)
    est <- s$estimates
    dg <- s$diagnostics
    dg$t_total <- eng$t_total
    draws <- stack_draws(eng$chains)
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      notes <- c(notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                dg$rhat_max, dg$ess_min))
    chains <- eng$chains
  } else {
    ml <- ms_fit_ml(des, settings)
    V <- ml$vcov
    ok_h <- !is.null(V) && all(is.finite(V)) && all(diag(V) > 0)
    boundary <- any(abs(ml$par) > 8)
    if (!ok_h) notes <- c(notes, "Hessian not invertible: some parameters are not identified.")
    if (boundary) notes <- c(notes, "Some link-scale estimates are beyond +/-8 (probabilities at 0 or 1): boundary estimate, do not interpret those.")
    pt <- matrix(ml$nimble_par, 1, dimnames = list(NULL, names(ml$nimble_par)))
    if (ok_h) {
      e <- eigen((V + t(V)) / 2, symmetric = TRUE)
      Z <- matrix(stats::rnorm(settings$n_draws * length(ml$par)), settings$n_draws)
      th <- sweep(Z %*% t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(ml$par))), 2, ml$par, "+")
      draws <- t(apply(th, 1, ms_theta_to_nimble, des = des))
    } else draws <- pt
    D <- ms_derive(draws, des); P0 <- ms_derive(pt, des)
    est <- data.frame(param = colnames(D), estimate = as.vector(P0), sd = apply(D, 2, stats::sd),
                      lower = apply(D, 2, stats::quantile, 0.025, na.rm = TRUE),
                      upper = apply(D, 2, stats::quantile, 0.975, na.rm = TRUE), row.names = NULL)
    dg <- list(converged = ml$convergence == 0 && ok_h && !boundary, rhat_max = NA_real_, ess_min = NA_real_,
               loglik = -ml$value)
    chains <- NULL
  }
  gof <- NULL
  if (isTRUE(settings$ppc_draws > 0)) {
    gof <- tryCatch(ms_ppc(draws, des, settings$ppc_draws), error = function(e) { notes <<- c(notes, paste("PPC failed:", conditionMessage(e))); NULL })
    if (!is.null(gof)) notes <- c(notes, ms_ppc_note(gof))
  }
  ucare <- NULL
  if (isTRUE(settings$ucare)) {
    ucare <- tryCatch(ms_ucare(des), error = function(e) { notes <<- c(notes, paste("U-CARE tests failed:", conditionMessage(e))); NULL })
    if (!is.null(ucare)) notes <- c(notes, ms_ucare_note(ucare, des))
  }
  dg$notes <- c(dg$notes, notes)
  list(estimates = est, diagnostics = dg,
       object = list(engine = engine, des = des, draws = draws, chains = chains, gof = gof, ucare = ucare,
                     settings = settings))
}

module <- ecomodule(
  id = "multistate",
  title = "Multistate and multievent capture-recapture",
  description = paste(
    "How likely are marked animals to survive in each state (site, breeding status, infection status), how often do they move",
    "or switch between states, and what is associated with those transitions? Arnason-Schwarz multistate model (Lebreton et al. 2009)",
    "for capture histories with the state recorded at each capture; multievent model (Pradel 2005) when the state of some detected",
    "animals could not be determined. Estimates survival per state, transition probabilities with covariate effects, state-specific",
    "detection, the probability that a detected animal's state is recorded, and the stationary state distribution."),
  defaults = list(
    params = list(survival_A = 0.8, survival_B = 0.65,
                  transition_A_to_B = 0.3, transition_B_to_A = 0.4,
                  detection_A = 0.6, detection_B = 0.5,
                  assignment_A = 0.7, assignment_B = 0.8, initial_A = 0.6, initial_B = 0.4,
                  survival_beta_size = 0.5,
                  transition_beta_temp_A_to_B = 0.8, transition_beta_temp_B_to_A = -0.5,
                  detection_beta_effort = 0.5),
    design = list(states = c("A", "B"), n_occasions = 8, n_marked = 40, model = "multievent",
                  survival_covs = "size", transition_covs = "temp", detection_covs = "effort")),
  simulate = ms_simulate, fit = ms_fit, settings = ms_settings, fast_settings = ms_fast,
  read_data = ms_read_data, check_data = ms_check_data, write_data = ms_write_data,
  interpret = ms_interpret, plots = ms_plots, figures = ms_figures, methods_text = ms_methods_text,
  power_params = c("survival_beta_size", "transition_beta_temp_A_to_B", "transition_beta_temp_B_to_A",
                   "detection_beta_effort"),
  engine = "nimble + nimbleEcology (dDHMMo); 'ml' = hand-written forward algorithm + optim",
  references = c(
    "Arnason AN (1973) The estimation of population size, migration rates and survival in a stratified population. Researches on Population Ecology 15:1-8. doi:10.1007/BF02510705",
    "Schwarz CJ, Schweigert JF, Arnason AN (1993) Estimating migration rates using tag-recovery data. Biometrics 49:177-193. doi:10.2307/2532612",
    "Brownie C, Hines JE, Nichols JD, Pollock KH, Hestbeck JB (1993) Capture-recapture studies for multiple strata including non-Markovian transitions. Biometrics 49:1173-1187. doi:10.2307/2532259",
    "Lebreton J-D, Nichols JD, Barker RJ, Pradel R, Spendelow JA (2009) Modeling individual animal histories with multistate capture-recapture models. Advances in Ecological Research 41:87-173. doi:10.1016/S0065-2504(09)00403-6",
    "Pradel R (2005) Multievent: an extension of multistate capture-recapture models to uncertain states. Biometrics 61:442-447. doi:10.1111/j.1541-0420.2005.00318.x",
    "Conn PB, Cooch EG (2009) Multistate capture-recapture analysis under imperfect state observation: an application to disease models. Journal of Applied Ecology 46:486-492. doi:10.1111/j.1365-2664.2008.01597.x",
    "Kendall WL (2009) One size does not fit all: adapting mark-recapture and occupancy models for state uncertainty. In: Thomson DL, Cooch EG, Conroy MJ (eds) Modeling Demographic Processes in Marked Populations. Springer, pp 765-780.",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.",
    "Zucchini W, MacDonald IL, Langrock R (2016) Hidden Markov Models for Time Series: An Introduction Using R, 2nd ed. CRC Press.",
    "Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS: A Hierarchical Perspective. Academic Press. Chapter 9: Estimation of survival and movement from capture-recapture data using multistate models.",
    "Turek D, de Valpine P, Paciorek CJ (2016) Efficient Markov chain Monte Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics 23:549-564. doi:10.1007/s10651-016-0353-z",
    "Goldstein B, Turek D, Ponisio L, de Valpine P (2024) nimbleEcology: Distributions for Ecological Models in nimble. R package version 0.5.0.",
    "Pradel R, Wintrebert CMA, Gimenez O (2003) A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model. Biometrics 59:43-53. doi:10.1111/1541-0420.00006",
    "Gimenez O, Lebreton J-D, Choquet R, Pradel R (2018) R2ucare: An R package to perform goodness-of-fit tests for capture-recapture models. Methods in Ecology and Evolution 9:1749-1754. doi:10.1111/2041-210X.13014",
    "Laake JL, Johnson DS, Conn PB (2013) marked: an R package for maximum likelihood and Markov Chain Monte Carlo analysis of capture-recapture data. Methods in Ecology and Evolution 4:885-890. doi:10.1111/2041-210X.12065"))
