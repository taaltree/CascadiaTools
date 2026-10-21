# fit(): resolve the engine(s), run them and assemble one estimates table.
#
# engine:
#   auto    (default) repeated visits -> gds; density covariates -> ds+hds;
#           otherwise ds
#   ds      conventional / multiple-covariate distance sampling (Distance::ds)
#   hds     hierarchical distance sampling (unmarked::distsamp) only
#   ds+hds  ds for density, abundance and detection; hds for covariate effects
#           on density (density_at_average_site, density_beta_*, rate_ratio_*)
#           and an HDS density over the study area (density_per_km2_hds)
#   gds     unmarked::gdistsamp (repeated visits, availability)
#   rlik    hand-written likelihoods (cross-check)

dist_settings <- list(engine = "auto", keys = c("hn", "hr"), adjustment = "cos", max_adjustments = 2,
                      hds_bins = 10, n_draws = 2000, gof_nboot = 200, K = NULL, rlik_key = "hn", model_average = TRUE,
                      variance = "analytic", nboot = 200)
dist_fast <- list(engine = "auto", keys = c("hn", "hr"), adjustment = "cos", max_adjustments = 2,
                  hds_bins = 10, n_draws = 500, gof_nboot = 0, K = NULL, rlik_key = "hn", model_average = TRUE,
                  variance = "analytic", nboot = 100)

dist_resolve_engine <- function(data, settings) {
  e <- settings$engine %||% "auto"
  if (e != "auto") return(e)
  if (data$n_visits > 1) return("gds")
  if (length(ds_terms(data$formulas$density))) return("ds+hds")
  "ds"
}

dist_fit <- function(data, settings) {
  settings <- merge_list(dist_settings, settings)
  data <- dist_prep(data)
  eng <- dist_resolve_engine(data, settings)
  parts <- list()
  if (eng %in% c("ds", "ds+hds")) parts$ds <- ds_engine_distance(data, settings)
  if (eng %in% c("hds", "ds+hds")) parts$hds <- ds_engine_unmarked(data, settings, gds = FALSE)
  if (eng == "gds") parts$gds <- ds_engine_unmarked(data, settings, gds = TRUE)
  if (eng == "rlik") parts$rlik <- ds_engine_rlik(data, settings)
  if (!length(parts)) stop("unknown engine '", eng, "'")
  if (eng == "ds+hds") {
    h <- parts$hds$estimates
    keep <- grepl("^(density_at_average_site|group_density_at_average_site|density_beta_|rate_ratio_)", h$param)
    hx <- h[h$param %in% c("density_per_km2", "abundance"), ]
    hx$param <- paste0(hx$param, "_hds")
    est <- rbind(parts$ds$estimates, h[keep, ], hx)
  } else est <- parts[[1]]$estimates
  rownames(est) <- NULL
  conv <- all(vapply(parts, function(p) isTRUE(p$converged), logical(1)))
  notes <- unlist(lapply(parts, function(p) p$notes))
  main <- parts[[1]]
  if (!is.null(main$gof) && !is.null(main$gof$p) && is.finite(main$gof$p) && main$gof$p < 0.05)
    notes <- c(notes, sprintf("Goodness-of-fit test rejects the detection function (p = %.3f): inspect the Q-Q plot and histogram for heaping or evasive movement.", main$gof$p))
  list(estimates = est,
       diagnostics = list(converged = conv, notes = notes, engine = eng),
       object = list(parts = parts, engine = eng, settings = settings),
       prepared = data[c("obs_w", "tz", "scaling", "st_terms", "dt_terms")])
}
