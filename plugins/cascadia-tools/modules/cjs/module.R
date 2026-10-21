# cjs: Cormack-Jolly-Seber models of apparent survival (phi) and recapture
# probability (p) from individual capture histories.
#   Engines: marked (maximum likelihood, default), nimble (Bayesian; marginal
#   likelihood via nimbleEcology::dCJS_vv, or the latent alive-state
#   formulation with settings formulation: latent).
#   Effects: time, linear trend, group (sex / age class), individual and
#   occasion covariates, transients (tsm), age at marking (age_class), trap
#   dependence (td), random time effects on survival (Nimble).
#   Goodness of fit: R2ucare (TEST2 / TEST3 components, c-hat).

for (f in c("prep.R", "engines.R", "fit.R", "sim.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

module <- ecomodule(
  id = "cjs",
  title = "Cormack-Jolly-Seber apparent survival",
  description = paste(
    "What proportion of marked animals survive (and stay in the study area) from one occasion to the next,",
    "does this apparent survival differ between groups (sex, age class) or change with individual traits (body mass) and",
    "environmental conditions (winter severity), and how likely is a live animal to be recaptured?",
    "Uses individual capture-recapture histories (Cormack-Jolly-Seber model), conditioning on first capture, so",
    "it estimates survival but not abundance. Apparent survival = true survival x site fidelity."),
  defaults = list(
    params = list(survival = 0.6, detection_prob = 0.5,
                  phi_beta_sex_male = 0.4, phi_beta_winter_severity = -0.5, phi_beta_body_mass = 0.3,
                  phi_beta_age_class_juvenile = -1, p_beta_sex_male = 0, p_beta_td = 0,
                  prop_transient = 0, temporal_sd = 0, p_temporal_sd = 0),
    design = list(n_occasions = 8, n_marked_per_occasion = 50, prop_male = 0.5, prop_juvenile = 0,
                  first_year = 2015, mass_mean = 21, mass_sd = 1.8, winter_mean = 45, winter_sd = 15,
                  phi_formula = "~ sex + winter_severity + body_mass", p_formula = "~ 1")),
  simulate = cjs_simulate, fit = cjs_fit, settings = cjs_settings, fast_settings = cjs_fast,
  read_data = cjs_read_data, check_data = cjs_check_data, write_data = cjs_write_data,
  interpret = cjs_interpret, plots = cjs_plots, figures = cjs_figures, methods_text = cjs_methods_text,
  power_params = c("phi_beta_sex_male", "phi_beta_winter_severity", "phi_beta_body_mass", "phi_beta_age_class_juvenile",
                   "p_beta_sex_male", "p_beta_td", "phi_beta_Time"),
  engine = "marked (ML), nimble + nimbleEcology (Bayesian, marginal or latent-state), R2ucare (GOF)",
  references = c(
    "Cormack RM (1964) Estimates of survival from the sighting of marked animals. Biometrika 51:429-438. doi:10.1093/biomet/51.3-4.429",
    "Jolly GM (1965) Explicit estimates from capture-recapture data with both death and immigration-stochastic model. Biometrika 52:225-247. doi:10.1093/biomet/52.1-2.225",
    "Seber GAF (1965) A note on the multiple-recapture census. Biometrika 52:249-259. doi:10.1093/biomet/52.1-2.249",
    "Lebreton J-D, Burnham KP, Clobert J, Anderson DR (1992) Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. Ecological Monographs 62:67-118. doi:10.2307/2937171",
    "Pradel R (1993) Flexibility in survival analysis from recapture data: handling trap-dependence. In: Lebreton J-D, North PM (eds) Marked Individuals in the Study of Bird Population, pp. 29-37. Birkhauser, Basel.",
    "Pradel R, Hines JE, Lebreton J-D, Nichols JD (1997) Capture-recapture survival models taking account of transients. Biometrics 53:60-72. doi:10.2307/2533097",
    "Pradel R, Gimenez O, Lebreton J-D (2005) Principles and interest of GOF tests for multistate capture-recapture models. Animal Biodiversity and Conservation 28:189-204.",
    "Gimenez O, Lebreton J-D, Choquet R, Pradel R (2018) R2ucare: An R package to perform goodness-of-fit tests for capture-recapture models. Methods in Ecology and Evolution 9:1749-1754. doi:10.1111/2041-210X.13014",
    "Choquet R, Lebreton J-D, Gimenez O, Reboulet A-M, Pradel R (2009) U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data. Ecography 32:1071-1074. doi:10.1111/j.1600-0587.2009.05968.x",
    "Laake JL, Johnson DS, Conn PB (2013) marked: an R package for maximum likelihood and Markov Chain Monte Carlo analysis of capture-recapture data. Methods in Ecology and Evolution 4:885-890. doi:10.1111/2041-210X.12065",
    "Gimenez O, Rossi V, Choquet R, Dehais C, Doris B, Varella H, Vila J-P, Pradel R (2007) State-space modelling of data on marked individuals. Ecological Modelling 206:431-438. doi:10.1016/j.ecolmodel.2007.03.040",
    "Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS: A Hierarchical Perspective, chapter 7 (Estimation of survival from capture-recapture data using the Cormack-Jolly-Seber model). Academic Press.",
    "Schaub M, Kery M (2022) Integrated Population Models: Theory and Ecological Applications with R and JAGS, chapter 3. Academic Press.",
    "Turek D, de Valpine P, Paciorek CJ (2016) Efficient Markov chain Monte Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics 23:549-564. doi:10.1007/s10651-016-0353-z",
    "Goldstein BR, Turek D, Ponisio L, de Valpine P (2024) nimbleEcology: Distributions for ecological models in nimble. R package version 0.5.0; https://CRAN.R-project.org/package=nimbleEcology",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487",
    "Burnham KP, Anderson DR (2002) Model Selection and Multimodel Inference, 2nd ed. Springer (quasi-likelihood c-hat adjustment)."))
