# ipm: integrated population model (Besbeas et al. 2002; Schaub & Abadi
# 2011; Schaub & Kery 2022, Integrated Population Models, chapters 5, 6, 9
# and 11). Female-based two-age-class pre-breeding-census model fitted
# jointly to population counts (state-space), capture-recapture data
# (age-dependent CJS m-arrays) and productivity data (Poisson), with
# immigration estimated as a latent rate. Engine: Nimble. Independent route:
# the same model in JAGS (book code) and the CJS component in `marked`.

for (f in c("model.R", "simulate.R", "derived.R", "io.R", "fit.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

ipm_defaults <- list(
  params = list(juvenile_survival = 0.25, adult_survival = 0.48,
                fecundity_yearling = 2.8, fecundity_adult = 3.6,
                immigration_rate = 0.12, recapture_prob = 0.6,
                juvenile_survival_sd = 0.3, adult_survival_sd = 0.2,
                fecundity_yearling_sd = 0.15, fecundity_adult_sd = 0.15,
                immigration_sd = 0.3, juvenile_survival_beta = -0.4),
  design = list(n_years = 15, N_initial = 100, marked_juveniles = 50, marked_adults = 20,
                prop_broods_monitored = 0.5, count_model = "poisson", count_sd = 10,
                immigration = TRUE, age_specific_fecundity = TRUE, temporal_variation = TRUE,
                covariate = TRUE, covariate_name = "winter_severity", sex_ratio = 0.5,
                first_year = 2001))

module <- ecomodule(
  id = "ipm",
  title = "Integrated population model (counts + capture-recapture + productivity)",
  description = paste(
    "How many breeding females are there each year, which vital rates (juvenile survival, adult survival,",
    "fecundity, immigration) drive the changes in population size, and how much does immigration contribute?",
    "An integrated population model fits one female-based, two-age-class population model jointly to annual",
    "counts, capture-recapture data of marked individuals and brood productivity data. This gives annual",
    "population size, growth rates, survival by age, fecundity and an immigration rate for which there are no",
    "direct data, and a life-table response experiment that apportions variation in population growth among",
    "the vital rates."),
  defaults = ipm_defaults,
  simulate = function(params, design, seed) { s <- ipm_simulate(params, design, seed); s[c("data", "truth")] },
  fit = ipm_fit, settings = ipm_settings, fast_settings = ipm_fast,
  read_data = ipm_read_data, check_data = ipm_check_data, write_data = ipm_write_data,
  interpret = ipm_interpret, figures = ipm_figures,
  plots = function(fit, data) ipm_figures(fit, data, fig_opts()),
  methods_text = ipm_methods_text,
  null_values = c(juvenile_survival_beta_winter_severity = 0),
  power_params = c("juvenile_survival_beta_winter_severity"),
  engine = "nimble (MCMC); cross-checks: JAGS (book code), marked (CJS maximum likelihood)",
  references = ipm_references)
