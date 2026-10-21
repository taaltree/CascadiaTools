# Graph utilities shared by both engines and the simulator: parse component
# models, build the edge list, enumerate directed paths, and compute direct /
# indirect / total effects (with Monte Carlo, delta-method or bootstrap CIs).

SEM_FAMILIES <- c("gaussian", "poisson", "nbinom", "binomial")

sem_canon_family <- function(f) {
  f <- tolower(f %||% "gaussian")
  f <- switch(f, normal = "gaussian", gauss = "gaussian", negbin = , negative.binomial = , nb = , nbinom2 = ,
              "negative binomial" = "nbinom", bernoulli = , logistic = "binomial", count = "poisson", f)
  if (!f %in% SEM_FAMILIES)
    stop(sprintf("family '%s' is not supported; use one of %s", f, paste(SEM_FAMILIES, collapse = ", ")))
  f
}

sem_default_link <- function(fam) switch(fam, gaussian = "identity", poisson = "log", nbinom = "log", binomial = "logit")

# models: list of list(formula = "y ~ a + b + (1|site)", family = "poisson", link = NULL)
sem_parse_models <- function(models) {
  if (is.null(models) || !length(models)) stop("no component models given (config `model:` list)")
  if (!is.null(models$formula)) models <- list(models)
  lapply(models, function(m) {
    if (is.character(m)) m <- list(formula = m)
    f <- stats::as.formula(m$formula)
    if (length(f) != 3) stop("component formula '", m$formula, "' needs a response: y ~ x")
    resp <- all.vars(f[[2]])
    if (length(resp) != 1 || !identical(deparse(f[[2]]), resp))
      stop(sprintf("response in '%s' must be a plain column name (make transformed columns with data: transform:)", m$formula))
    bars <- lme4::findbars(f)
    fixed <- lme4::nobars(f)
    tl <- attr(stats::terms(fixed), "term.labels")
    if (any(grepl(":", tl))) stop(sprintf("interactions are not supported in component '%s' (indirect effects through interactions are not defined); see README Pitfalls", m$formula))
    bad <- tl[!vapply(tl, function(t) identical(make.names(t), t), logical(1))]
    if (length(bad)) stop(sprintf("predictor '%s' in '%s' must be a plain column name; create transformed columns with data: transform:", bad[1], m$formula))
    groups <- unique(unlist(lapply(bars, function(b) all.vars(b[[3]]))))
    if (length(bars) && !all(vapply(bars, function(b) identical(deparse(b[[2]]), "1"), logical(1))))
      stop(sprintf("only random intercepts (1|group) are supported in '%s'", m$formula))
    fam <- sem_canon_family(m$family)
    link <- m$link %||% sem_default_link(fam)
    list(formula = f, text = paste(deparse(f), collapse = ""), response = resp, preds = tl, groups = groups,
         family = fam, link = link)
  })
}

# Response first, in topological order (errors on cycles)
sem_topo_order <- function(comps) {
  resp <- vapply(comps, `[[`, "", "response")
  if (anyDuplicated(resp)) stop("response '", resp[duplicated(resp)][1], "' has more than one component model; combine them into one formula")
  done <- character(); ord <- integer()
  repeat {
    ready <- which(!seq_along(comps) %in% ord & vapply(comps, function(c) all(setdiff(c$preds, done) %in% setdiff(c$preds, resp)), logical(1)))
    if (!length(ready)) break
    ord <- c(ord, ready); done <- c(done, resp[ready])
  }
  if (length(ord) < length(comps)) stop("the model has a feedback loop (cycle) among: ", paste(setdiff(resp, done), collapse = ", "), ". Piecewise SEM needs an acyclic model.")
  comps[ord]
}

# All directed paths from `from` to `to` in an edge list (data.frame from, to)
sem_paths <- function(edges, from, to, max_paths = 500) {
  out <- list()
  walk <- function(node, path) {
    if (length(out) >= max_paths) return()
    if (node == to) { out[[length(out) + 1]] <<- path; return() }
    nx <- edges$to[edges$from == node]
    for (k in nx) if (!k %in% path) walk(k, c(path, k))
  }
  walk(from, from)
  out
}

# Effects for pairs. edges: data.frame(from, to, key) where key indexes theta.
# pairs: data.frame(from, to). Returns a function theta -> named vector of
# direct, indirect and total effects (sum over paths of products).
sem_effect_fun <- function(edges, pairs, prefix = "") {
  specs <- lapply(seq_len(nrow(pairs)), function(i) {
    ps <- sem_paths(edges, pairs$from[i], pairs$to[i])
    keys <- lapply(ps, function(p) vapply(seq_len(length(p) - 1), function(j)
      edges$key[edges$from == p[j] & edges$to == p[j + 1]][1], ""))
    list(from = pairs$from[i], to = pairs$to[i], paths = ps, keys = keys,
         direct = vapply(keys, length, 1L) == 1)
  })
  nm <- unlist(lapply(specs, function(s) {
    b <- paste0(s$from, "_to_", s$to)
    paste0(prefix, c("direct_", "indirect_", "total_"), b)
  }))
  f <- function(theta) {
    v <- unlist(lapply(specs, function(s) {
      pr <- vapply(s$keys, function(k) prod(theta[k]), 1)
      d <- sum(pr[s$direct]); ind <- sum(pr[!s$direct])
      c(d, ind, d + ind)
    }))
    stats::setNames(v, nm)
  }
  attr(f, "specs") <- specs
  f
}

# Default effects of interest: every (exogenous or endogenous) source -> any
# endogenous sink reachable by at least one indirect path.
sem_default_pairs <- function(edges) {
  nodes <- unique(c(edges$from, edges$to))
  pr <- list()
  for (a in nodes) for (b in setdiff(nodes, a)) {
    ps <- sem_paths(edges, a, b)
    if (any(vapply(ps, length, 1L) > 2)) pr[[length(pr) + 1]] <- data.frame(from = a, to = b)
  }
  if (!length(pr)) return(data.frame(from = character(), to = character()))
  do.call(rbind, pr)
}

# Intermediate nodes on all indirect paths for a pair
sem_mediators <- function(spec) unique(unlist(lapply(spec$paths, function(p) if (length(p) > 2) p[-c(1, length(p))])))

# CI rows for a vector-valued function g(theta) using draws / delta method.
# draws: matrix (ndraw x length(theta)) with colnames = names(theta)
sem_ci_rows <- function(g, theta, method = "montecarlo", vcov = NULL, draws = NULL, prob = 0.95) {
  est <- g(theta)
  a <- (1 - prob) / 2
  if (method == "delta") {
    k <- length(theta); J <- matrix(0, length(est), k)
    for (j in seq_len(k)) {
      h <- 1e-6 * max(1, abs(theta[j])); tp <- theta; tm <- theta
      tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
      J[, j] <- (g(tp) - g(tm)) / (2 * h)
    }
    se <- sqrt(pmax(0, diag(J %*% vcov %*% t(J))))
    z <- stats::qnorm(1 - a)
    return(data.frame(param = names(est), estimate = unname(est), sd = se, lower = est - z * se, upper = est + z * se, row.names = NULL))
  }
  vals <- do.call(rbind, lapply(seq_len(nrow(draws)), function(i) g(draws[i, ])))
  data.frame(param = names(est), estimate = unname(est), sd = apply(vals, 2, stats::sd),
             lower = apply(vals, 2, stats::quantile, a, na.rm = TRUE),
             upper = apply(vals, 2, stats::quantile, 1 - a, na.rm = TRUE), row.names = NULL)
}

# Multivariate normal draws without MASS dependency issues
sem_rmvn <- function(n, mu, S) {
  if (!length(mu)) return(matrix(0, n, 0))
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(mu))
  z <- matrix(stats::rnorm(n * length(mu)), n)
  sweep(z %*% t(L), 2, mu, "+")
}

fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
fmt_sig <- function(x, d = 3) trimws(formatC(signif(x, d), digits = d, format = "fg"))
