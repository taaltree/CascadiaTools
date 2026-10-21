# ecotools module `sem`: structural equation models for ecology.
#
# Data structure (what simulate()$data and read_data() both return):
#   list(df   = data.frame, one row per sampling unit (plot, site, animal),
#        spec = list(engine = "piecewise" | "lavaan",
#                    models = list(list(formula, family, link)),   # piecewise
#                    lavaan = "<lavaan syntax>", group = <column>,  # lavaan
#                    effects = list(list(from, to)),  # indirect/total effects to report
#                    labels, units, estimator, drop_missing))
# Engines: piecewiseSEM (component (G)LM(M)s + d-separation; default) and
# lavaan (latent variables, multigroup, correlated errors). See README.md.

for (f in c("graph.R", "simulate.R", "engine_piecewise.R", "engine_lavaan.R", "output.R"))
  sys.source(file.path(MODULE_DIR, f), envir = environment())

SEM_SETTINGS <- list(ci_method = "montecarlo", n_mc = 4000, n_boot = 500, seed = 1,
                     dsep_direction = "causal", dsep = TRUE, r2 = TRUE, estimator = "ML")

sem_fit <- function(data, settings) {
  eco_require(c("piecewiseSEM", "lavaan", "lme4", "MASS"))
  st <- merge_list(SEM_SETTINGS, settings)
  eng <- data$spec$engine %||% "piecewise"
  if (eng == "lavaan") sem_fit_lavaan(data, st) else sem_fit_piecewise(data, st)
}

module <- ecomodule(
  id = "sem",
  title = "Structural equation models (piecewise SEM and lavaan)",
  description = paste(
    "How do several ecological variables influence one another as a network: which links are supported by the data,",
    "how large is each direct path, and how much of a driver's association with an outcome runs indirectly through",
    "mediators? Fits piecewise SEMs (one GLM/GLMM per response, tested with Shipley's d-separation test) or",
    "lavaan SEMs with latent variables and multiple groups, and reports path coefficients (raw and standardised),",
    "indirect and total effects with confidence intervals, and global goodness of fit."),
  defaults = list(
    params = list(fire_to_canopy = -0.6, canopy_to_shrub = -0.5, moisture_to_shrub = 0.4,
                  shrub_to_birds = 0.4, fire_to_birds = -0.25, int_canopy = 0, int_shrub = 0,
                  int_birds = 1.5, sd_canopy = 1, sd_shrub = 1),
    design = list(
      n = 100, engine = "piecewise",
      model = list(list(formula = "canopy ~ fire", family = "gaussian"),
                   list(formula = "shrub ~ canopy + moisture", family = "gaussian"),
                   list(formula = "birds ~ shrub + fire", family = "poisson")),
      effects = list(list(from = "fire", to = "birds"), list(from = "moisture", to = "birds"),
                     list(from = "fire", to = "shrub")))),
  simulate = sem_simulate, fit = sem_fit,
  settings = SEM_SETTINGS,
  fast_settings = merge_list(SEM_SETTINGS, list(n_mc = 2000, r2 = FALSE)),
  read_data = sem_read_data, check_data = sem_check_data,
  interpret = sem_interpret,
  figures = sem_figures,
  plots = function(fit, data) {
    o <- fig_opts()
    figs <- sem_figures(fit, data, o)
    figs <- figs[setdiff(names(figs), "effects_panel")]
    lapply(stats::setNames(names(figs), names(figs)), function(n) if (inherits(figs[[n]], "ggplot")) apply_fig_opts(figs[[n]], n, merge_list(o, list(base_size = 11))) else figs[[n]])
  },
  methods_text = sem_methods,
  power_params = c("path_fire_to_birds", "path_shrub_to_birds", "path_fire_to_canopy", "path_moisture_to_shrub",
                   "indirect_fire_to_birds", "indirect_moisture_to_birds", "total_fire_to_birds"),
  engine = "piecewiseSEM (component lm/glm/lme4/MASS models + d-separation) or lavaan (ML/MLR)",
  references = c(
    "Grace JB (2006). Structural Equation Modeling and Natural Systems. Cambridge University Press.",
    "Grace JB, Anderson TM, Olff H, Scheiner SM (2010). On the specification of structural equation models for ecological systems. Ecological Monographs 80(1):67-87.",
    "Shipley B (2000). A new inferential test for path models based on directed acyclic graphs. Structural Equation Modeling 7(2):206-218.",
    "Shipley B (2009). Confirmatory path analysis in a generalized multilevel context. Ecology 90(2):363-368.",
    "Shipley B, Douma JC (2020). Generalized AIC and chi-squared statistics for path models consistent with directed acyclic graphs. Ecology 101(3):e02960.",
    "Lefcheck JS (2016). piecewiseSEM: Piecewise structural equation modelling in R for ecology, evolution, and systematics. Methods in Ecology and Evolution 7(5):573-579.",
    "Lefcheck JS. Structural Equation Modeling in R for Ecology and Evolution (online book, chapters 'Local estimation' and 'Coefficients'). https://jslefche.github.io/sem_book/",
    "Grace JB, Johnson DJ, Lefcheck JS, Byrnes JEK (2018). Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.",
    "Rosseel Y (2012). lavaan: An R package for structural equation modeling. Journal of Statistical Software 48(2):1-36. doi:10.18637/jss.v048.i02",
    "Preacher KJ, Selig JP (2012). Advantages of Monte Carlo confidence intervals for indirect effects. Communication Methods and Measures 6(2):77-98. doi:10.1080/19312458.2012.679848",
    "Hu L, Bentler PM (1999). Cutoff criteria for fit indexes in covariance structure analysis: conventional criteria versus new alternatives. Structural Equation Modeling 6(1):1-55.",
    "Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software 67(1):1-48."))
