# occu-spatial: spatial and multi-species occupancy via spOccupancy
# (Doser et al. 2022). Variants chosen with settings$model:
#   spPGOcc   single species, NNGP exponential spatial random effect (default for single)
#   PGOcc     single species, non-spatial
#   msPGOcc   multi-species community model, non-spatial
#   sfMsPGOcc multi-species spatial factor model (default for multi)
#   intPGOcc  integrated occupancy, several data sources (default for integrated)
#   spIntPGOcc  integrated + spatial random effect (experimental; not validated)
# The data type (single / multi / integrated) comes from the data
# (species column -> multi, source column -> integrated).

for (f in c("sim.R", "fit.R", "io.R", "report.R")) {
  ff <- file.path(MODULE_DIR, "R", f)
  if (file.exists(ff)) sys.source(ff, envir = environment())
}

module <- ecomodule(
  id = "occu-spatial",
  title = "Spatial, multi-species and integrated occupancy (spOccupancy)",
  description = paste(
    "Where does a species (or a community of species) occur, which habitat covariates are associated with occurrence,",
    "and how much residual spatial clustering is there, accounting for imperfect detection? Variants: single-species",
    "spatial occupancy with a nearest-neighbour Gaussian process (spPGOcc), multi-species community occupancy (msPGOcc),",
    "multi-species spatial factor occupancy (sfMsPGOcc) and integrated occupancy combining several survey methods (intPGOcc)."),
  defaults = list(
    params = list(occupancy_prob = 0.5, psi_beta_forest = 0.8,
                  detection_prob = 0.4, p_beta_effort = 0.5,
                  spatial_variance = 1.5, effective_range_km = 6),
    design = list(type = "single", n_sites = 200, n_surveys = 4, extent_km = 20,
                  prop_missing = 0.1, cov_spatial_frac = 0.5, cov_range_km = 8)),
  simulate = simulate_occu_spatial,
  fit = fit_occu_spatial,
  settings = list(model = "auto", n_iter = 50000, n_burn = 20000, n_thin = 20, n_chains = 3,
                  batch_length = 25, n_neighbors = 15, n_factors = 1, gof = TRUE,
                  converge_on = "all", seed = 1),
  fast_settings = list(model = "auto", n_iter = 15000, n_burn = 5000, n_thin = 10, n_chains = 3,
                       batch_length = 25, n_neighbors = 10, n_factors = 1, gof = FALSE,
                       converge_on = "core", seed = 1),
  read_data = if (exists("read_occu_spatial")) read_occu_spatial else NULL,
  check_data = if (exists("check_occu_spatial")) check_occu_spatial else NULL,
  interpret = if (exists("interpret_occu_spatial")) interpret_occu_spatial else NULL,
  figures = if (exists("figures_occu_spatial")) figures_occu_spatial else NULL,
  plots = if (exists("plots_occu_spatial")) plots_occu_spatial else NULL,
  methods_text = if (exists("methods_occu_spatial")) methods_occu_spatial else NULL,
  power_params = c("psi_beta_forest", "p_beta_effort", "comm_psi_beta_forest"),
  engine = "spOccupancy (Polya-Gamma Gibbs samplers, NNGP spatial random effects)",
  references = c(
    "Doser JW, Finley AO, Kery M, Zipkin EF (2022) spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. Methods in Ecology and Evolution 13:1670-1678. doi:10.1111/2041-210X.13897",
    "Doser JW, Finley AO, Banerjee S (2023) Joint species distribution models with imperfect detection for high-dimensional spatial data. Ecology 104(9):e4137.",
    "Doser JW, Leuenberger W, Sillett TS, Hallworth MT, Zipkin EF (2022) Integrated community occupancy models: a framework to assess occurrence and biodiversity dynamics using multiple data sources. Methods in Ecology and Evolution 13. doi:10.1111/2041-210X.13811",
    "Datta A, Banerjee S, Finley AO, Gelfand AE (2016) Hierarchical nearest-neighbor Gaussian process models for large geostatistical datasets. Journal of the American Statistical Association 111:800-812.",
    "Polson NG, Scott JG, Windle J (2013) Bayesian inference for logistic models using Polya-Gamma latent variables. Journal of the American Statistical Association 108:1339-1349.",
    "Dorazio RM, Royle JA (2005) Estimating size and composition of biological communities by modeling the occurrence of species. Journal of the American Statistical Association 100:389-398.",
    "MacKenzie DI, Nichols JD, Lachman GB, Droege S, Royle JA, Langtimm CA (2002) Estimating site occupancy rates when detection probabilities are less than one. Ecology 83:2248-2255.",
    "Banerjee S, Carlin BP, Gelfand AE (2014) Hierarchical Modeling and Analysis for Spatial Data, 2nd ed. CRC Press.")
)
