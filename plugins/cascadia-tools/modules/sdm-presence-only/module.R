# sdm-presence-only: species distribution models from presence-only records,
# done as inhomogeneous Poisson point processes (Warton & Shepherd 2010;
# Renner et al. 2015), with model-based correction for observer bias
# (Warton et al. 2013; Fithian et al. 2015) and optional integration of
# survey data with non-detections (Koshkina et al. 2017; Dorazio 2014).
#
# Engines (settings$engine):
#   ppm         down-weighted Poisson regression (glm) on grid quadrature (default without surveys)
#   integrated  joint likelihood: presence-only IPP + survey occupancy (cloglog, detection p)
#   maxnet      Maxent via glmnet (Phillips et al. 2017), percentile-bootstrap intervals
#   survey_only the survey part alone (for comparison)

for (f in c("simulate.R", "prepare.R", "engines.R", "fit.R", "post.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

module <- ecomodule(
  id = "sdm-presence-only",
  title = "Species distribution from presence-only records (point process, Maxent, integrated)",
  description = paste(
    "Where in the landscape is the species most likely to occur, and which habitat covariates is that associated with,",
    "using presence-only records (museum/herbarium specimens, citizen-science sightings, incidental observations)?",
    "Records are modelled as an inhomogeneous Poisson point process whose intensity depends on habitat, thinned by",
    "observer effort that depends on a sampling-bias covariate (e.g. distance to road). Presence-only data give a",
    "RELATIVE intensity map and habitat effects; adding survey sites with non-detections (integrated model) also",
    "gives absolute density, expected abundance, occupancy probability and detection probability."),
  defaults = list(
    params = list(density_at_mean_covariates = 0.5,
                  intensity_beta_forest_cover_pct = 0.8, intensity_beta_elevation_m = -0.5,
                  bias_beta_dist_road_km = -1.0, detection_prob = 0.5),
    design = list(nx = 60, ny = 60, cell_km = 1, n_presence = 200,
                  survey_sites = 0, survey_visits = 3, survey_area_km2 = 1,
                  bias_habitat_cor = 0.6, range_km = 6, n_roads = 4)),
  simulate = sdm_simulate, fit = sdm_fit,
  settings = sdm_settings, fast_settings = sdm_fast,
  read_data = sdm_read_data, check_data = sdm_check_data, write_data = sdm_write_data,
  interpret = sdm_interpret,
  plots = function(fit, data) sdm_figures(fit, data, fig_opts()),
  figures = sdm_figures, methods_text = sdm_methods_text,
  power_params = c("intensity_beta_forest_cover_pct", "intensity_beta_elevation_m", "bias_beta_dist_road_km"),
  engine = "glm (down-weighted Poisson regression), optim (integrated likelihood), maxnet",
  references = sdm_references)
