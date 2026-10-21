# TEMPLATE module: logistic regression of presence on one covariate.
# Copy this folder to modules/<your-id>/ and replace each piece. It shows the
# full contract (see core/R/contract.R) with the smallest possible model.
#
# Data structure (what simulate()$data and read_data() both return):
#   list(df = data.frame(y = 0/1, x = numeric))

simulate <- function(params, design, seed) {
  set.seed(seed)
  x <- stats::rnorm(design$n_sites)
  p <- stats::plogis(params$beta0 + params$beta1 * x)
  y <- stats::rbinom(design$n_sites, 1, p)
  list(data = list(df = data.frame(y = y, x = x)),
       truth = c(beta0 = params$beta0, beta1 = params$beta1))
}

fit <- function(data, settings) {
  m <- stats::glm(y ~ x, family = stats::binomial, data = data$df)
  co <- summary(m)$coefficients
  est <- wald_rows(co[, 1], co[, 2], c("beta0", "beta1"))
  list(estimates = est,
       diagnostics = list(converged = m$converged, notes = character()),
       object = m)
}

read_data <- function(config) {
  df <- utils::read.csv(file.path(config$.dir, config$data$file))
  list(df = data.frame(y = df[[config$data$response]], x = df[[config$data$covariate]]))
}

check_data <- function(data) {
  pr <- character()
  if (!all(data$df$y %in% c(0, 1))) pr <- c(pr, "response must be 0/1")
  if (anyNA(data$df)) pr <- c(pr, "missing values present")
  pr
}

interpret <- function(fit, data) {
  b <- fit$estimates[fit$estimates$param == "beta1", ]
  sprintf("Each one-unit increase in the covariate multiplies the odds of presence by %.2f (95%% CI %.2f to %.2f).",
          exp(b$estimate), exp(b$lower), exp(b$upper))
}

# Publication figure suite (docs/FIGURES.md). opts = fig_opts(config$figures):
# use fig_label() so users can rename anything, and draw covariates in
# original units.
figures <- function(fit, data, opts = fig_opts()) {
  m <- fit$object
  co <- fit$estimates
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  xs <- seq(min(data$df$x), max(data$df$x), length.out = 100)
  pr <- stats::predict(m, newdata = data.frame(x = xs), se.fit = TRUE)
  z <- stats::qnorm(1 - (1 - ci) / 2)
  pred <- data.frame(x = xs, estimate = stats::plogis(pr$fit),
                     lower = stats::plogis(pr$fit - z * pr$se.fit), upper = stats::plogis(pr$fit + z * pr$se.fit))
  list(
    effect_x = with_caption(
      fig_cond_effect(pred, fig_label("x", "Covariate", opts), fig_label("y", "Probability of presence", opts),
                      obs = if (isTRUE(opts$show_data)) data.frame(x = data$df$x, y = data$df$y) else NULL, ylim01 = TRUE),
      sprintf("Predicted probability of presence (line) with %d%% confidence band; points are observed sites.", round(100 * ci))),
    coefficients = with_caption(
      fig_coef(co[co$param != "beta0", ], null_value = 1, transform = exp, log_x = TRUE,
               xlab = "Odds ratio per unit covariate (95% CI)", opts = opts),
      "Odds ratio for the covariate; dashed line = no effect.")
  )
}

module <- ecomodule(
  id = "_template",
  title = "Template: logistic regression",
  description = "Does the probability that a species is present at a site change with a site covariate?",
  defaults = list(params = list(beta0 = 0, beta1 = 0.5), design = list(n_sites = 200)),
  simulate = simulate, fit = fit, read_data = read_data, check_data = check_data,
  interpret = interpret, figures = figures, power_params = "beta1", engine = "stats::glm",
  references = "McCullagh P, Nelder JA (1989) Generalized Linear Models. Chapman & Hall."
)
