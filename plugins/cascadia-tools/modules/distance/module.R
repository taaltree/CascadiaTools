# distance: distance sampling for density and abundance.
#   ds      conventional (CDS) and multiple-covariate (MCDS) distance sampling with
#           Distance::ds: half-normal / hazard-rate keys, cosine or Hermite
#           adjustments chosen by AIC, line or point transects, exact or binned
#           distances, truncation, strata, clusters; density and abundance with
#           the Fewster et al. (2009) encounter-rate variance; Cramer-von Mises GOF.
#   hds     hierarchical distance sampling (Royle et al. 2004; unmarked::distsamp):
#           site covariates on density.
#   gds     unmarked::gdistsamp (Chandler et al. 2011): repeated visits, availability.
#   rlik    hand-written likelihoods (independent cross-check).

for (f in c("detfun.R", "prep.R", "sim.R", "engine_ds.R", "engine_hds.R", "engine_rlik.R",
            "fit.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

module <- ecomodule(
  id = "distance",
  title = "Distance sampling (conventional, multiple-covariate and hierarchical)",
  description = paste(
    "How many animals are there per square kilometre and in the whole study area, and is density associated with habitat,",
    "when observers walk lines or stand at points and record the distance to each animal (or group) they see, missing more of them farther away?",
    "Conventional and multiple-covariate distance sampling (Buckland et al. 2001, 2015) fit a detection function to the distances",
    "and give density, abundance and the probability of detecting an animal inside the truncation distance;",
    "hierarchical distance sampling (Royle et al. 2004) adds transect-level covariates on density;",
    "repeated visits (Chandler et al. 2011) also estimate availability."),
  defaults = list(
    params = list(group_density_at_average_site = 4, density_beta_habitat = 0.4, sigma_m = 150, hr_shape = 3,
                  detection_key = "hn", mean_group_size = 1, availability_prob = 0.6, detection_beta_observer = 0),
    design = list(survey = "line", n_transects = 30, transect_length_km = 3, cell_width_km = 2, point_cell_km2 = 0.25,
                  truncation_m = 400, distance_type = "exact", n_bins = 8, n_visits = 1, n_strata = 1, observers = FALSE,
                  responsive_move_m = 0, heap_zero_m = 0, round_m = 0, fixed_positions = FALSE)),
  simulate = dist_simulate, fit = dist_fit, settings = dist_settings, fast_settings = dist_fast,
  read_data = dist_read_data, check_data = dist_check_data,
  interpret = dist_interpret, plots = dist_plots, figures = dist_figures, methods_text = dist_methods_text,
  power_params = c("density_beta_habitat", "rate_ratio_habitat", "detection_beta_observerB"),
  null_values = c(density_beta_habitat = 0, rate_ratio_habitat = 1, detection_beta_observerB = 0),
  engine = "Distance::ds (CDS/MCDS), unmarked::distsamp / gdistsamp (HDS), hand-written likelihood (cross-check)",
  references = dist_references)
