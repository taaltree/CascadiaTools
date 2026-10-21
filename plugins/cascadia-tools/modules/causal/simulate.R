# Structural causal model (SCM) simulator for the causal module.
#
# Every node except the exposure and the outcome is linear-Gaussian in its
# parents and standardised (mean 0, SD 1 in the population), so the path
# coefficients b_<from>__<to> are standardised path coefficients:
#     v = sum_p b_pv * p + e_v,  Var(e_v) = 1 - Var(sum_p b_pv * p)
# Exposure: continuous (as above) or binary, X ~ Bernoulli(plogis(a_x + lin)).
# Outcome: link(E[Y]) = eta = intercept + sum_p b_py * p, with the family's
#   noise (gaussian sigma, poisson, negbin theta, binomial trials, Gamma shape).
# Children of the outcome use the standardised outcome (Y - E[Y]) / SD[Y].
# Panel designs: nodes in design[["unit_level"]] are drawn once per unit (e.g. an
#   unmeasured site-level confounder U), nodes in design[["time_level"]] once per
#   time period; everything else once per row.
# Optional natural units: design[["node_scale"]] = list(node = c(mean, sd)).
#
# Truth (cached per params/design): the causal estimands are computed from
# the SCM by intervention, with a large Monte Carlo sample (4e5) using
# common random numbers, so they are exact up to ~0.2% Monte Carlo error:
#   beta_<x>   total effect on the link scale per unit of x = sum over
#              directed paths of the products of path coefficients (tau).
#              Given only for identity and log links (a logit coefficient
#              is non-collapsible and has no single population value).
#   ate_<x>    E[Y(x + delta)] - E[Y(x)] (continuous: every unit's exposure
#              shifted by `contrast`) or E[Y(1)] - E[Y(0)] (binary).
#   <ratio>_<x>  E[Y(x + delta)] / E[Y(x)] (rate/risk/mean ratio, non-identity links)
#   mean_<y>_if_<x>_0/1  potential-outcome means (binary exposure)
#   direct_*   the same for the controlled direct effect (edge x -> y only).

.causal_cache <- new.env()

causal_family_info <- function(family, link = NULL) {
  family <- tolower(family %||% "gaussian")
  if (family %in% c("nb", "nbinom", "nbinom2", "negative_binomial", "negative.binomial")) family <- "negbin"
  if (family == "gamma") family <- "Gamma"
  def <- switch(family, gaussian = "identity", poisson = "log", negbin = "log", binomial = "logit", Gamma = "log",
                stop("family must be one of gaussian, poisson, negbin, binomial, Gamma (got '", family, "')"))
  link <- link %||% def
  ratio <- switch(family, gaussian = NA_character_, poisson = "rate_ratio", negbin = "rate_ratio",
                  binomial = "risk_ratio", Gamma = "mean_ratio")
  if (link == "identity") ratio <- if (family == "gaussian") NA_character_ else ratio
  linkinv <- switch(link, identity = identity, log = exp, logit = stats::plogis, probit = stats::pnorm,
                    stop("link '", link, "' not supported"))
  list(family = family, link = link, ratio = ratio, linkinv = linkinv)
}

# Names of the estimands (shared by simulate() truth and fit()).
causal_estimand_names <- function(x, y, fam, exposure_type, prefix = "") {
  xs <- cn_clean(x); ys <- cn_clean(y)
  n <- c(beta = paste0(prefix, "beta_", xs), ate = paste0(prefix, "ate_", xs))
  if (!is.na(fam$ratio)) n["ratio"] <- paste0(prefix, fam$ratio, "_", xs)
  if (exposure_type == "binary" && prefix == "") {
    n["mean0"] <- paste0("mean_", ys, "_if_", xs, "_0"); n["mean1"] <- paste0("mean_", ys, "_if_", xs, "_1")
  }
  n
}

causal_coef_matrix <- function(g, params) {
  e <- causal_edges(g)
  nodes <- names(g)
  B <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  miss <- character()
  for (k in seq_len(nrow(e))) {
    nm <- paste0("b_", e$from[k], "__", e$to[k])
    v <- params[[nm]]
    if (is.null(v)) { v <- params$default_path; if (is.null(v)) miss <- c(miss, nm) }
    if (!is.null(v)) B[e$from[k], e$to[k]] <- as.numeric(v)
  }
  if (length(miss)) stop("simulate(): missing path coefficients (add to params): ", paste(miss, collapse = ", "),
                         ". Or give params$default_path.")
  B
}

causal_topo <- function(g) {
  o <- unlist(dagitty::topologicalOrdering(g))
  names(sort(o))
}

# Generate the SCM. `noise` = list of standard normal vectors per node (and
# `u_y` uniform for the outcome); `calib` = intercepts/noise SDs; `do` =
# list(node = value vector) intervention; returns all nodes + eta.
causal_generate <- function(g, B, x, y, fam, params, design, calib, n, noise, do = NULL) {
  ord <- causal_topo(g)
  v <- list()
  xb <- identical(design[["exposure_type"]], "binary")
  for (node in ord) {
    pa <- dagitty::parents(g, node)
    lin <- if (length(pa)) Reduce(`+`, lapply(pa, function(p) B[p, node] * (if (p == y) v$.y_std else v[[p]]))) else rep(0, n)
    if (!is.null(do[[node]])) { v[[node]] <- rep_len(do[[node]], n); next }
    if (node == y) {
      eta <- (params$intercept %||% 0) + lin
      v$.eta <- eta
      mu <- fam$linkinv(eta)
      tr <- design[["trials"]] %||% 1
      yy <- switch(fam$family,
        gaussian = mu + (params$sigma %||% 1) * noise[[node]],
        poisson = stats::qpois(noise$.u_y, mu),
        negbin = stats::qnbinom(noise$.u_y, size = params$theta %||% 2, mu = mu),
        binomial = stats::qbinom(noise$.u_y, tr, mu),
        Gamma = stats::qgamma(noise$.u_y, shape = params$shape %||% 4, rate = (params$shape %||% 4) / mu))
      v[[node]] <- yy
      v$.y_std <- (yy - calib$y_mean) / calib$y_sd
    } else if (node == x && xb) {
      p <- stats::plogis((params$exposure_intercept %||% 0) + lin - calib$center[[node]])
      v[[node]] <- as.numeric(noise$.u_x < p)
    } else {
      v[[node]] <- lin - calib$center[[node]] + calib$sd[[node]] * noise[[node]]
    }
  }
  v
}

causal_noise <- function(g, n, design, seed, n_units = NULL, unit = NULL, n_times = NULL, time = NULL) {
  set.seed(seed)
  nodes <- names(g)
  nz <- lapply(stats::setNames(nodes, nodes), function(v) {
    if (!is.null(unit) && v %in% unlist(design[["unit_level"]])) stats::rnorm(n_units)[unit]
    else if (!is.null(time) && v %in% unlist(design[["time_level"]])) stats::rnorm(n_times)[time]
    else stats::rnorm(n)
  })
  nz$.u_x <- stats::runif(n); nz$.u_y <- stats::runif(n)
  nz
}

# Calibration (centering, noise SDs, outcome moments) and truth, cached.
causal_truth <- function(params, design) {
  key <- paste(deparse(list(params, design[setdiff(names(design), c("n_sites", "n_units", "n_times"))])), collapse = "")
  if (!is.null(.causal_cache[[key]])) return(.causal_cache[[key]])
  eco_require("dagitty")
  g <- causal_parse(causal_dag_string(design[["dag"]], design[["edges"]]))
  x <- design[["exposure"]]; y <- design[["outcome"]]
  fam <- causal_family_info(design[["family"]], design[["link"]])
  B <- causal_coef_matrix(g, params)
  ord <- causal_topo(g)
  N <- 400000
  nz <- causal_noise(g, N, design, seed = 987654)
  # sequential calibration: standardise each Gaussian node in topological order
  calib <- list(center = list(), sd = list(), y_mean = 0, y_sd = 1)
  v <- list()
  xb <- identical(design[["exposure_type"]], "binary")
  for (node in ord) {
    pa <- dagitty::parents(g, node)
    lin <- if (length(pa)) Reduce(`+`, lapply(pa, function(p) B[p, node] * (if (p == y) v$.y_std else v[[p]]))) else rep(0, N)
    if (node == y) {
      eta <- (params$intercept %||% 0) + lin
      mu <- fam$linkinv(eta)
      tr <- design[["trials"]] %||% 1
      yy <- switch(fam$family,
        gaussian = mu + (params$sigma %||% 1) * nz[[node]],
        poisson = stats::qpois(nz$.u_y, mu),
        negbin = stats::qnbinom(nz$.u_y, size = params$theta %||% 2, mu = mu),
        binomial = stats::qbinom(nz$.u_y, tr, mu),
        Gamma = stats::qgamma(nz$.u_y, shape = params$shape %||% 4, rate = (params$shape %||% 4) / mu))
      calib$y_mean <- mean(yy); calib$y_sd <- stats::sd(yy); v$.eta <- eta
      v[[node]] <- yy; v$.y_std <- (yy - calib$y_mean) / calib$y_sd
    } else if (node == x && xb) {
      calib$center[[node]] <- mean(lin); calib$sd[[node]] <- 0
      v[[node]] <- as.numeric(nz$.u_x < stats::plogis((params$exposure_intercept %||% 0) + lin - mean(lin)))
    } else {
      vl <- stats::var(lin)
      if (vl > 0.98) stop(sprintf("simulate(): path coefficients into '%s' explain %.0f%% of its variance (must be < 98%%); make them smaller.", node, 100 * vl))
      calib$center[[node]] <- mean(lin); calib$sd[[node]] <- sqrt(1 - vl)
      v[[node]] <- lin - mean(lin) + calib$sd[[node]] * nz[[node]]
    }
  }
  # total effect on the link scale (sum of path products), via (I - B)^-1
  Tm <- solve(diag(nrow(B)) - B)
  tau <- Tm[x, y]; bdir <- B[x, y]
  sc <- design[["node_scale"]][[x]]
  sx <- if (!xb && !is.null(sc)) sc[2] else 1
  delta <- if (xb) 1 else (design[["contrast"]] %||% (if (!is.null(sc)) sc[2] else 1)) / sx  # in SD units
  eta <- v$.eta
  ex_type <- if (xb) "binary" else "continuous"
  tr <- list()
  mk <- function(shift, prefix, b) {
    nm <- causal_estimand_names(x, y, fam, ex_type, prefix)
    out <- c()
    if (fam$link %in% c("identity", "log")) out[nm["beta"]] <- b / sx
    if (xb) { hi <- eta + b * (1 - v[[x]]); lo <- eta + b * (0 - v[[x]]) } else { hi <- eta + b * delta; lo <- eta }
    mh <- mean(fam$linkinv(hi)); ml <- mean(fam$linkinv(lo))
    out[nm["ate"]] <- mh - ml
    if (!is.na(fam$ratio)) out[nm["ratio"]] <- mh / ml
    if (xb && prefix == "") { out[nm["mean0"]] <- ml; out[nm["mean1"]] <- mh }
    out
  }
  truth <- mk(delta, "", tau)
  eff <- design[["effect"]] %||% "total"
  if (eff %in% c("direct", "both")) truth <- c(truth, mk(delta, "direct_", bdir))
  res <- list(truth = truth, calib = calib, B = B, g = g, fam = fam, tau = tau, direct = bdir, delta_sd = delta)
  .causal_cache[[key]] <- res
  res
}

causal_simulate <- function(params, design, seed) {
  eco_require("dagitty")
  tt <- causal_truth(params, design)
  g <- tt$g; x <- design[["exposure"]]; y <- design[["outcome"]]
  panel <- !is.null(design[["n_units"]])
  if (panel) {
    nu <- design[["n_units"]]; nt <- design[["n_times"]] %||% 5
    n <- nu * nt; unit <- rep(seq_len(nu), each = nt); time <- rep(seq_len(nt), nu)
  } else { n <- design[["n_sites"]] %||% 150; unit <- NULL; time <- NULL; nu <- NULL; nt <- NULL }
  nz <- causal_noise(g, n, design, seed, nu, unit, nt, time)
  v <- causal_generate(g, tt$B, x, y, tt$fam, params, design, tt$calib, n, nz)
  lat <- union(dagitty::latents(g), unlist(design[["latent"]]))
  keep <- setdiff(names(g), lat)
  df <- as.data.frame(v[keep])
  for (nd in setdiff(keep, y)) {
    sc <- design[["node_scale"]][[nd]]
    if (!is.null(sc) && !(nd == x && identical(design[["exposure_type"]], "binary"))) df[[nd]] <- sc[1] + sc[2] * df[[nd]]
  }
  if (!is.null(design[["round"]])) for (nd in names(design[["round"]])) if (nd %in% names(df)) df[[nd]] <- round(df[[nd]], design[["round"]][[nd]])
  if (panel) { df <- cbind(data.frame(unit = sprintf("U%03d", unit), time = time), df) }
  tr <- design[["trials"]] %||% 1
  if (tt$fam$family == "binomial" && tr > 1) df$trials <- tr
  spec <- list(dag = causal_dag_string(design[["dag"]], design[["edges"]]), exposure = x, outcome = y,
               family = tt$fam$family, link = design[["link"]], effect = design[["effect"]] %||% "total",
               contrast = if (identical(design[["exposure_type"]], "binary")) NULL else (design[["contrast"]] %||% (design[["node_scale"]][[x]][2] %||% 1)),
               adjust = design[["adjust"]], panel = design[["panel"]] %||% "none",
               unit = if (panel) "unit", time = if (panel) "time",
               unit_constant = design[["unit_constant"]], time_constant = design[["time_constant"]],
               trials = if (tt$fam$family == "binomial" && tr > 1) "trials")
  truth <- tt$truth
  if (isTRUE(design[["whatif_truth"]])) {
    # truths for the what-if rows: every strategy targets the same causal
    # effect, so a biased strategy must fail the recovery gate.
    id <- causal_identify(spec, names(df))
    st <- causal_strategies(id$g, x, y, id$total$set, id$observed, spec[["panel"]])
    tags <- setdiff(st$tag, "dag")
    if (spec[["panel"]] != "none") tags <- c(tags, "no_fixed_effects")
    nm <- causal_estimand_names(x, y, tt$fam, if (identical(design[["exposure_type"]], "binary")) "binary" else "continuous")
    truth <- c(truth, stats::setNames(rep(truth[[nm[["ate"]]]], length(tags)), paste0("whatif_", tags, "__", nm[["ate"]])))
  }
  list(data = list(df = df, spec = spec), truth = truth)
}
