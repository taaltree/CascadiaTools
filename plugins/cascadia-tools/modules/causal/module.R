# ecotools module `causal`: causal effects from observational data with DAGs.
#
# Data structure (what simulate()$data and read_data() both return):
#   list(df   = data.frame, one row per site (or per site x time for panels);
#               DAG node names are column names,
#        spec = list(dag, exposure, outcome, family, link, trials, offset,
#                    offset_per, contrast, exposure_reference, effect,
#                    adjust, adjust_direct, panel, unit, time,
#                    unit_constant, time_constant, labels, units,
#                    response_label))
# See README.md for every key.
#
# Pipeline: DAG -> roles and adjustment sets (dagitty) -> GLM with the
# adjustment set (optionally unit / unit+time fixed effects) -> average
# causal effect by g-computation (marginaleffects) -> what-if comparison
# under wrong adjustment sets -> local tests of the DAG, DHARMa, E-value.

for (f in c("dag.R", "simulate.R", "engine.R", "output.R")) sys.source(file.path(MODULE_DIR, f), envir = environment())

# Default structural causal model: does riparian canopy cover increase
# juvenile trout abundance in streams? (standardised path coefficients)
CAUSAL_STREAM_DAG <- "dag {
  elevation -> canopy
  elevation -> water_temp
  agriculture -> canopy
  agriculture -> sediment
  sediment -> trout
  width -> canopy
  width -> trout
  canopy -> water_temp
  water_temp -> trout
  canopy -> trout
  canopy -> algae
  trout -> algae
}"

causal_defaults <- list(
  params = list(
    b_elevation__canopy = 0.4, b_elevation__water_temp = -0.6,
    b_agriculture__canopy = -0.5, b_agriculture__sediment = 0.6,
    b_sediment__trout = -0.4, b_width__canopy = -0.3, b_width__trout = 0.3,
    b_canopy__water_temp = -0.4, b_water_temp__trout = -0.4,
    b_canopy__trout = 0.2, b_canopy__algae = -0.4, b_trout__algae = 0.5,
    intercept = log(5), theta = 3, sigma = 1, exposure_intercept = 0),
  design = list(dag = CAUSAL_STREAM_DAG, exposure = "canopy", outcome = "trout",
                family = "negbin", n_sites = 120, contrast = 1, effect = "total",
                exposure_type = "continuous"))

module <- ecomodule(
  id = "causal",
  title = "Causal effects from observational data (DAG + g-computation)",
  description = paste(
    "What is the causal effect of an exposure (habitat, management, a stressor) on an ecological outcome, estimated from",
    "observational data? You draw your causal assumptions as a DAG; the module finds which variables must (and must not) be",
    "adjusted for, checks the DAG against the data, estimates the average causal effect on the response scale with a GLM of",
    "any family (or a fixed-effects panel estimator when sites are repeatedly measured and have unmeasured site-level",
    "confounders), and shows how wrong adjustment (omitting confounders, adjusting for mediators or colliders) would change the answer."),
  defaults = causal_defaults,
  simulate = causal_simulate, fit = causal_fit,
  settings = CAUSAL_DEFAULT_SETTINGS,
  fast_settings = merge_list(CAUSAL_DEFAULT_SETTINGS, list(whatif = FALSE, local_tests = FALSE, diagnostics = FALSE)),
  read_data = causal_read_data, check_data = causal_check_data,
  interpret = causal_interpret, figures = causal_figures, plots = causal_plots,
  methods_text = causal_methods_text,
  power_params = c("beta_canopy", "ate_canopy"),
  null_values = c(beta_canopy = 0, ate_canopy = 0, rate_ratio_canopy = 1),
  engine = "dagitty (identification) + stats::glm / MASS::glm.nb + marginaleffects (g-computation); sandwich (cluster-robust SE)",
  references = c(
    "Pearl J (2009) Causality: Models, Reasoning, and Inference, 2nd edn. Cambridge University Press.",
    "Arif S, MacNeil MA (2023) Applying the structural causal model framework for observational causal inference in ecology. Ecological Monographs 93(1): e1554. doi:10.1002/ecm.1554",
    "Laubach ZM, Murray EJ, Hoke KL, Safran RJ, Perng W (2021) A biologist's guide to model selection and causal inference. Proceedings of the Royal Society B 288: 20202815. doi:10.1098/rspb.2020.2815",
    "Byrnes JEK, Dee LE (2025) Causal inference with observational data and unobserved confounding variables. Ecology Letters 28(1): e70023. doi:10.1111/ele.70023",
    "Textor J, van der Zander B, Gilthorpe MS, Liskiewicz M, Ellison GTH (2016) Robust causal inference using directed acyclic graphs: the R package 'dagitty'. International Journal of Epidemiology 45(6): 1887-1894. doi:10.1093/ije/dyw341",
    "Cinelli C, Forney A, Pearl J (2022) A crash course in good and bad controls. Sociological Methods & Research 53(3): 1071-1104. doi:10.1177/00491241221099552",
    "Henckel L, Perkovic E, Maathuis MH (2022) Graphical criteria for efficient total effect estimation via adjustment in causal linear models. Journal of the Royal Statistical Society B 84: 579-599. doi:10.1111/rssb.12451",
    "Arel-Bundock V, Greifer N, Heiss A (2024) How to interpret statistical models using marginaleffects for R and Python. Journal of Statistical Software 111(9): 1-32. doi:10.18637/jss.v111.i09",
    "Shipley B (2000) A new inferential test for path models based on directed acyclic graphs. Structural Equation Modeling 7(2): 206-218.",
    "VanderWeele TJ, Ding P (2017) Sensitivity analysis in observational research: introducing the E-value. Annals of Internal Medicine 167: 268-274. doi:10.7326/M16-2607",
    "Mundlak Y (1978) On the pooling of time series and cross section data. Econometrica 46: 69-85. doi:10.2307/1913646",
    "Hartig F (2022) DHARMa: residual diagnostics for hierarchical (multi-level/mixed) regression models. R package.")
)
