# Detection functions and their integrals (all distances in metres).
#   half-normal  g(x) = exp(-x^2 / (2 sigma^2))
#   hazard-rate  g(x) = 1 - exp(-(x / sigma)^(-b))          (Hayes & Buckland 1983)
# Line transects: P_a = (1/w) int_0^w g(x) dx; ESW = w * P_a.
# Point transects: P_a = (2/w^2) int_0^w r g(r) dr; EDR = w * sqrt(P_a).

ds_g <- function(x, key, sigma, b = NULL) {
  if (key == "hn") exp(-x^2 / (2 * sigma^2))
  else if (key == "hr") 1 - exp(-(pmax(x, 1e-12) / sigma)^(-b))
  else if (key == "unif") rep(1, length(x))
  else stop("unknown detection key ", key)
}

# Average detection probability within truncation w, vectorised over sigma
# (and b). Half-normal uses the closed form; hazard-rate a 401-point Simpson
# rule (relative error < 1e-6 for the shapes seen in practice).
ds_pa <- function(key, sigma, w, survey = "line", b = NULL, left = 0) {
  n <- max(length(sigma), length(b))
  sigma <- rep_len(sigma, n)
  if (!is.null(b)) b <- rep_len(b, n)
  if (key == "hn" && left == 0) {
    if (survey == "line") return(sigma * sqrt(2 * pi) * (stats::pnorm(w / sigma) - 0.5) / w)
    return(2 * sigma^2 * (1 - exp(-w^2 / (2 * sigma^2))) / w^2)
  }
  if (key == "unif") return(rep(1, n))
  m <- 400
  x <- seq(left, w, length.out = m + 1)
  wt <- c(1, rep(c(4, 2), length.out = m - 1), 1) * (w - left) / (3 * m)
  out <- vapply(seq_len(n), function(i) {
    g <- ds_g(x, key, sigma[i], b[i])
    if (survey == "line") sum(wt * g) / (w - left) else sum(wt * 2 * x * g) / (w^2 - left^2)
  }, 1)
  out
}

# Probability mass of each distance bin among detections (cutpoints in m),
# relative to the covered area: line: int_bin g / w, point: int_bin 2 r g / w^2.
ds_bin_probs <- function(key, sigma, cut, survey = "line", b = NULL) {
  w <- max(cut)
  vapply(seq_len(length(cut) - 1), function(j) {
    a <- cut[j]; z <- cut[j + 1]
    if (key == "hn" && survey == "line")
      return(sigma * sqrt(2 * pi) * (stats::pnorm(z / sigma) - stats::pnorm(a / sigma)) / w)
    if (key == "hn" && survey == "point")
      return(2 * sigma^2 * (exp(-a^2 / (2 * sigma^2)) - exp(-z^2 / (2 * sigma^2))) / w^2)
    m <- 60
    x <- seq(a, z, length.out = m + 1)
    wt <- c(1, rep(c(4, 2), length.out = m - 1), 1) * (z - a) / (3 * m)
    g <- ds_g(x, key, sigma, b)
    if (survey == "line") sum(wt * g) / w else sum(wt * 2 * x * g) / w^2
  }, 1)
}

ds_esw <- function(pa, w, survey) if (survey == "line") pa * w else w * sqrt(pa)
