# Simulate data from a recursive path model (a DAG of GLM/GLMM components)
# and compute the true value of every estimand the fit reports.
#
# design$model: generative components, list(list(formula, family, link)).
#   Exogenous nodes (never a response) are drawn N(mean_<x>, sd_<x>).
# params (flat, so `power --what params --vary x_to_y=...` works):
#   <x>_to_<y>  path coefficient (link scale of y)    int_<y>  intercept
#   sd_<y>      residual SD (gaussian)                theta_<y> NB size
#   mean_<x>, sd_<x>  exogenous mean / SD             re_sd_<y> random-intercept SD
# design$fit_model: components to FIT (piecewise) when they differ from the
#   generative model (negative controls); design$lavaan: lavaan syntax to fit
#   with engine = "lavaan"; design$latent: generative nodes not observed.
# design$n_groups: number of groups for (1|group) terms (plots spread evenly).

sem_par <- function(params, key, default) {
  v <- params[[key]]
  if (is.null(v)) default else as.numeric(v)
}

sem_generate <- function(params, comps, n, n_groups = NULL) {
  comps <- sem_topo_order(comps)
  resp <- vapply(comps, `[[`, "", "response")
  exo <- setdiff(unique(unlist(lapply(comps, `[[`, "preds"))), resp)
  groups <- unique(unlist(lapply(comps, `[[`, "groups")))
  df <- data.frame(row.names = seq_len(n))
  for (g in groups) {
    ng <- n_groups %||% max(2, round(n / 10))
    df[[g]] <- factor(sprintf("%s%03d", g, sort(rep(seq_len(ng), length.out = n))))
  }
  for (x in exo) df[[x]] <- stats::rnorm(n, sem_par(params, paste0("mean_", x), 0), sem_par(params, paste0("sd_", x), 1))
  etas <- list()
  for (c in comps) {
    y <- c$response
    eta <- rep(sem_par(params, paste0("int_", y), 0), n)
    for (x in c$preds) eta <- eta + sem_par(params, paste0(x, "_to_", y), 0) * df[[x]]
    for (g in c$groups) {
      u <- stats::rnorm(nlevels(df[[g]]), 0, sem_par(params, paste0("re_sd_", y), 0.5))
      eta <- eta + u[as.integer(df[[g]])]
    }
    etas[[y]] <- eta
    df[[y]] <- switch(c$family,
      gaussian = eta + stats::rnorm(n, 0, sem_par(params, paste0("sd_", y), 1)),
      poisson = stats::rpois(n, exp(eta)),
      nbinom = stats::rnbinom(n, size = sem_par(params, paste0("theta_", y), 2), mu = exp(eta)),
      binomial = stats::rbinom(n, 1, if (c$link == "probit") stats::pnorm(eta) else stats::plogis(eta)))
  }
  attr(df, "eta") <- etas
  df
}

# Population scaling factors for standardised coefficients, matching the
# piecewiseSEM conventions (Lefcheck 2016; Grace et al. 2018): gaussian sd(y);
# binomial latent-theoretic sqrt(var(eta) + pi^2/3) (logit) or + 1 (probit);
# poisson / NB observation-empirical sd(eta) / cor(y, mu).
sem_sd_y <- function(y, eta, family, link) {
  switch(family,
    gaussian = stats::sd(y),
    binomial = sqrt(stats::var(eta) + if (link == "probit") 1 else pi^2 / 3),
    poisson = , nbinom = stats::sd(eta) / stats::cor(y, exp(eta)))
}

.sem_truth_cache <- new.env(parent = emptyenv())

sem_population_sd <- function(params, comps, n, n_groups, n_big = 2e5) {
  key <- paste(deparse(list(params, lapply(comps, function(c) c[c("text", "family", "link")]), n, n_groups)), collapse = "")
  if (!is.null(.sem_truth_cache[[key]])) return(.sem_truth_cache[[key]])
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  set.seed(424242)
  ngb <- if (is.null(n_groups)) NULL else max(2, round(n_big * n_groups / n))
  big <- sem_generate(params, comps, n_big, ngb)
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  eta <- attr(big, "eta")
  sds <- vapply(names(big)[vapply(big, is.numeric, logical(1))], function(v) stats::sd(big[[v]]), 1)
  sdy <- list()
  for (c in comps) sdy[[c$response]] <- sem_sd_y(big[[c$response]], eta[[c$response]], c$family, c$link)
  out <- list(sd_x = sds, sd_y = unlist(sdy))
  .sem_truth_cache[[key]] <- out
  out
}

sem_edges_from_comps <- function(comps, params = NULL) {
  e <- do.call(rbind, lapply(comps, function(c) if (length(c$preds))
    data.frame(from = c$preds, to = c$response, stringsAsFactors = FALSE)))
  if (is.null(e)) e <- data.frame(from = character(), to = character())
  e$key <- paste0(e$from, "_to_", e$to)
  if (!is.null(params)) e$value <- vapply(e$key, function(k) sem_par(params, k, 0), 1)
  e
}

# Edge naming used in the estimates table: path_ for structural paths,
# loading_ for latent -> indicator (lavaan measurement model).
sem_edge_param <- function(from, to, latent = character(), indicators = character()) {
  ifelse(from %in% latent & to %in% indicators, paste0("loading_", from, "_", to), paste0("path_", from, "_to_", to))
}

# Number of sampling units. `n_units` is the key to use in YAML configs: YAML
# 1.1 reads a bare key `n` as the boolean FALSE. `n` (CLI --vary n=...) is
# used when n_units is absent.
sem_design_n <- function(design) {
  if (!is.null(design[["FALSE"]])) stop("design key `n` in a YAML config is read as FALSE (YAML 1.1); write n_units: instead")
  as.integer(design$n_units %||% design$n)
}

sem_simulate <- function(params, design, seed) {
  set.seed(seed)
  gen <- sem_parse_models(design$model)
  n <- sem_design_n(design)
  df <- sem_generate(params, gen, n, design$n_groups)
  latent <- as.character(unlist(design$latent))
  attr(df, "eta") <- NULL
  df <- df[, setdiff(names(df), latent), drop = FALSE]
  engine <- design$engine %||% "piecewise"
  spec <- list(engine = engine, effects = design$effects, labels = design$labels, units = design$units,
               group = design$group, estimator = design$estimator, simulated = TRUE)
  if (engine == "lavaan") {
    spec$lavaan <- design$lavaan
  } else {
    spec$models <- design$fit_model %||% design$model
  }
  data <- list(df = df, spec = spec)
  list(data = data, truth = sem_truth(params, design, gen))
}

sem_truth <- function(params, design, gen = sem_parse_models(design$model)) {
  latent <- as.character(unlist(design$latent))
  gen_e <- sem_edges_from_comps(gen, params)
  # fitted graph (edges only)
  fit_e <- if ((design$engine %||% "piecewise") == "lavaan") sem_lavaan_edges_syntax(design$lavaan)
           else sem_edges_from_comps(sem_parse_models(design$fit_model %||% design$model))
  indicators <- if (length(latent)) gen_e$to[gen_e$from %in% latent & !gen_e$to %in% unique(c(gen_e$from))] else character()
  ps <- sem_population_sd(params, gen, sem_design_n(design), design$n_groups)
  sdx <- function(v) ps$sd_x[[v]]
  fam <- stats::setNames(vapply(gen, `[[`, "", "family"), vapply(gen, `[[`, "", "response"))
  std_fac <- function(from, to) sdx(from) / ps$sd_y[[to]]
  all_e <- unique(rbind(gen_e[, c("from", "to")], fit_e[, c("from", "to")]))
  tr <- c()
  for (i in seq_len(nrow(all_e))) {
    f <- all_e$from[i]; t <- all_e$to[i]
    v <- sem_par(params, paste0(f, "_to_", t), 0)
    if (!any(gen_e$from == f & gen_e$to == t)) v <- 0
    if (!any(fit_e$from == f & fit_e$to == t)) next   # only edges that are estimated
    nm <- sem_edge_param(f, t, latent, indicators)
    tr[nm] <- v
    tr[paste0("std_", sub("^path_", "", nm))] <- v * std_fac(f, t)
  }
  pairs <- sem_pairs_from_spec(design$effects, if (is.null(fit_e$op)) fit_e else fit_e[fit_e$op == "~", ])
  if (nrow(pairs)) {
    gen_e$std <- gen_e$value * mapply(std_fac, gen_e$from, gen_e$to)
    th <- stats::setNames(gen_e$value, gen_e$key); ths <- stats::setNames(gen_e$std, gen_e$key)
    tr <- c(tr, sem_effect_fun(gen_e, pairs)(th), sem_effect_fun(gen_e, pairs, "std_")(ths))
    fe <- paste(fit_e$from, fit_e$to)
    for (i in seq_len(nrow(pairs))) if (!paste(pairs$from[i], pairs$to[i]) %in% fe)
      tr <- tr[!names(tr) %in% paste0(c("direct_", "std_direct_"), pairs$from[i], "_to_", pairs$to[i])]
  }
  tr
}

sem_pairs_from_spec <- function(effects, edges) {
  if (is.null(effects) || !length(effects)) return(sem_default_pairs(edges))
  if (is.data.frame(effects)) return(effects[, c("from", "to")])
  do.call(rbind, lapply(effects, function(e) data.frame(from = e$from, to = e$to, stringsAsFactors = FALSE)))
}
