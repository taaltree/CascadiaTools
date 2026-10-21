# DAG handling for the causal module: parsing, node roles, adjustment sets
# (dagitty), the optimal adjustment set (Henckel et al. 2022) and the
# "what if" adjustment strategies used for the educational comparison table.
#
# Conventions
#   - DAG node names ARE the data column names (after the optional
#     `columns:` map in config has been applied by read_data()).
#   - A node that is not a column of the data is unobserved (latent).
#   - Latent nodes listed in `unit_constant` / `time_constant` are absorbed
#     by unit / time fixed effects when a panel estimator is used, so for
#     identification they are treated as "adjusted for".

cn_clean <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Build the dagitty object from a string ("dag { a -> b ... }", or just the
# edge statements) or from an edge list (character "a -> b" or data.frame
# from/to).
causal_dag_string <- function(dag = NULL, edges = NULL) {
  if (!is.null(dag)) {
    s <- paste(dag, collapse = "\n")
    if (!grepl("^\\s*(dag|pdag|graph)\\s*\\{", s)) s <- paste0("dag {\n", s, "\n}")
    return(s)
  }
  if (is.null(edges)) stop("config needs `dag:` (dagitty syntax) or `edges:` (list of 'from -> to')")
  if (is.data.frame(edges)) e <- paste(edges$from, "->", edges$to)
  else e <- unlist(edges)
  paste0("dag {\n", paste(e, collapse = "\n"), "\n}")
}

causal_parse <- function(dag_string) {
  eco_require("dagitty")
  g <- tryCatch(dagitty::dagitty(dag_string), error = function(e)
    stop("could not read the DAG: ", conditionMessage(e), ". Use dagitty syntax, e.g. 'dag { forest -> birds; elevation -> forest; elevation -> birds }'"))
  if (!identical(dagitty::graphType(g), "dag")) stop("the graph must be a DAG (only -> arrows); got a ", dagitty::graphType(g))
  if (!dagitty::isAcyclic(g)) stop("the DAG has a cycle (a variable is its own cause). Feedback loops must be unrolled in time (x_t -> y_t+1).")
  g
}

# Set the latent flags: nodes absent from the data, plus declared latents,
# minus nodes absorbed by fixed effects (if `absorb` is TRUE).
causal_set_latent <- function(g, observed, absorbed = character(), absorb = TRUE) {
  nodes <- names(g)
  lat <- union(dagitty::latents(g), setdiff(nodes, observed))
  lat <- setdiff(lat, intersect(lat, observed))
  if (absorb) lat <- setdiff(lat, absorbed)
  dagitty::latents(g) <- intersect(lat, nodes)
  g
}

causal_edges <- function(g) {
  e <- dagitty::edges(g)
  data.frame(from = as.character(e$v), to = as.character(e$w), stringsAsFactors = FALSE)
}

# Roles of every node relative to exposure x and outcome y.
causal_roles <- function(g, x, y) {
  nodes <- names(g)
  dx <- setdiff(dagitty::descendants(g, x), x)
  ay <- setdiff(dagitty::ancestors(g, y), y)
  ax <- setdiff(dagitty::ancestors(g, x), x)
  dy <- setdiff(dagitty::descendants(g, y), y)
  role <- stats::setNames(rep("other", length(nodes)), nodes)
  role[intersect(ax, ay)] <- "confounder (cause of exposure and outcome)"
  role[setdiff(ax, ay)] <- "cause of exposure only"
  role[setdiff(intersect(setdiff(ay, ax), setdiff(nodes, dx)), x)] <- "cause of outcome only"
  role[intersect(dx, ay)] <- "mediator (on a causal path)"
  role[intersect(dx, dy)] <- "collider / descendant of outcome"
  role[setdiff(dy, dx)] <- "descendant of outcome"
  role[setdiff(intersect(dx, setdiff(nodes, ay)), dy)] <- "descendant of exposure only"
  role[x] <- "exposure"; role[y] <- "outcome"
  role
}

# Proper-causal-path nodes (excluding x) and the forbidden set.
causal_cn <- function(g, x, y) {
  intersect(setdiff(dagitty::descendants(g, x), x), dagitty::ancestors(g, y))
}

# Optimal adjustment set O(x, y) of Henckel, Perkovic & Maathuis (2022):
# parents of the causal nodes, minus the forbidden set. NULL if it contains
# an unobserved node or is not a valid adjustment set.
causal_optimal_set <- function(g, x, y) {
  cn <- causal_cn(g, x, y)
  forb <- union(unlist(lapply(cn, function(v) dagitty::descendants(g, v))), x)
  O <- setdiff(unique(unlist(dagitty::parents(g, cn))), forb)
  if (length(intersect(O, dagitty::latents(g)))) return(NULL)
  if (!isTRUE(dagitty::isAdjustmentSet(g, O, x, y))) return(NULL)
  sort(O)
}

causal_minimal_sets <- function(g, x, y, effect = c("total", "direct")) {
  effect <- match.arg(effect)
  s <- dagitty::adjustmentSets(g, x, y, type = "minimal", effect = effect)
  lapply(s, function(v) sort(as.character(v)))
}

# Pick one set: user list > "optimal" > smallest minimal set (ties broken
# towards sets whose members are causes of the outcome, which tend to give
# more precise estimates; Henckel et al. 2022).
causal_choose_set <- function(g, x, y, adjust = NULL, effect = "total") {
  sets <- causal_minimal_sets(g, x, y, effect)
  notes <- character()
  if (is.character(adjust) && !identical(adjust, "minimal") && !identical(adjust, "optimal")) {
    ok <- if (effect == "total") isTRUE(dagitty::isAdjustmentSet(g, adjust, x, y))
          else any(vapply(sets, function(s) setequal(s, adjust), logical(1)))
    return(list(set = sort(adjust), sets = sets, how = "user", valid = ok, notes = notes))
  }
  if (!length(sets)) return(list(set = NULL, sets = sets, how = "none", valid = FALSE, notes = notes))
  if (identical(adjust, "optimal") && effect == "total") {
    O <- causal_optimal_set(g, x, y)
    if (!is.null(O)) return(list(set = O, sets = sets, how = "optimal", valid = TRUE, notes = notes))
    notes <- c(notes, "The optimal adjustment set includes an unmeasured variable, so the smallest minimal set was used instead.")
  }
  pa_y <- dagitty::parents(g, y)
  size <- vapply(sets, length, 1L)
  near <- vapply(sets, function(s) sum(s %in% pa_y), 1L)
  k <- order(size, -near)[1]
  list(set = sets[[k]], sets = sets, how = "minimal", valid = TRUE, notes = notes)
}

# "What if" strategies: what a researcher might adjust for instead, and
# whether the DAG says each identifies the total effect.
causal_strategies <- function(g, x, y, chosen, observed, panel = "none", max_rows = 14) {
  roles <- causal_roles(g, x, y)
  obs <- setdiff(intersect(names(g), observed), c(x, y))
  lat <- dagitty::latents(g)
  valid <- function(s) isTRUE(dagitty::isAdjustmentSet(g, s, x, y))
  out <- list()
  add <- function(tag, label, set, kind, short = label) {
    set <- sort(unique(set))
    key <- paste(set, collapse = ",")
    if (any(vapply(out, function(o) identical(o$key, key), logical(1)))) return(invisible())
    out[[length(out) + 1]] <<- list(tag = tag, label = label, set = set, key = key, kind = kind,
                                    valid = valid(set), short = short)
  }
  add("dag", "DAG adjustment set (used for the main estimate)", chosen, "dag", "DAG set (main estimate)")
  add("crude", "No adjustment (crude association)", character(), "crude", "No adjustment (crude)")
  for (v in intersect(chosen, setdiff(chosen, lat))) {
    s <- setdiff(chosen, v)
    add(paste0("omit_", cn_clean(v)), sprintf("Omit %s (%s)", v, if (grepl("^confounder", roles[v])) "a confounder" else "it blocks a back-door path"), s, "omit",
        sprintf("Omit %s", v))
  }
  for (v in obs[grepl("^mediator", roles[obs])])
    add(paste0("add_mediator_", cn_clean(v)), sprintf("Also adjust for mediator %s (over-control)", v), c(chosen, v), "mediator", sprintf("+ mediator %s", v))
  for (v in obs[grepl("collider|descendant of outcome", roles[obs])])
    add(paste0("add_collider_", cn_clean(v)),
        sprintf("Also adjust for %s (%s)", v, if (grepl("collider", roles[v])) "collider: caused by exposure and outcome" else "a consequence of the outcome"),
        c(chosen, v), "collider", sprintf("+ %s %s", if (grepl("collider", roles[v])) "collider" else "outcome descendant", v))
  for (v in setdiff(obs, chosen)) {
    if (roles[v] %in% c("cause of exposure only")) add(paste0("add_instrument_", cn_clean(v)), sprintf("Also adjust for %s (cause of exposure only: valid but less precise)", v), c(chosen, v), "instrument", sprintf("+ exposure-only cause %s", v))
    if (roles[v] %in% c("cause of outcome only")) add(paste0("add_precision_", cn_clean(v)), sprintf("Also adjust for %s (cause of outcome only: valid, often more precise)", v), c(chosen, v), "precision", sprintf("+ outcome-only cause %s", v))
    if (roles[v] %in% c("confounder (cause of exposure and outcome)", "other") && !valid(c(chosen, v)))
      add(paste0("add_bias_", cn_clean(v)), sprintf("Also adjust for %s (opens a biasing path, e.g. M-bias)", v), c(chosen, v), "collider", sprintf("+ %s (opens a path)", v))
  }
  sets <- causal_minimal_sets(g, x, y)
  for (i in seq_along(sets)) add(paste0("alt_set_", i), sprintf("Alternative valid set {%s}", paste(sets[[i]], collapse = ", ")), sets[[i]], "alternative", sprintf("Valid set {%s}", paste(sets[[i]], collapse = ", ")))
  O <- causal_optimal_set(g, x, y)
  if (!is.null(O)) add("optimal", sprintf("Optimal set {%s} (Henckel et al. 2022)", paste(O, collapse = ", ")), O, "alternative", sprintf("Optimal set {%s}", paste(O, collapse = ", ")))
  add("all_measured", "Adjust for every measured variable (kitchen sink)", obs, "kitchen", "All measured variables")
  out <- out[seq_len(min(length(out), max_rows))]
  df <- data.frame(tag = vapply(out, `[[`, "", "tag"), label = vapply(out, `[[`, "", "label"),
                   kind = vapply(out, `[[`, "", "kind"), short = vapply(out, `[[`, "", "short"), valid = vapply(out, `[[`, TRUE, "valid"),
                   stringsAsFactors = FALSE)
  df$set <- lapply(out, `[[`, "set")
  df
}

# Full identification summary used by fit() and simulate().
causal_identify <- function(spec, observed) {
  g0 <- causal_parse(causal_dag_string(spec[["dag"]], spec[["edges"]]))
  x <- spec[["exposure"]]; y <- spec[["outcome"]]
  miss <- setdiff(c(x, y), names(g0))
  if (length(miss)) stop("exposure/outcome not in the DAG: ", paste(miss, collapse = ", "))
  panel <- spec[["panel"]] %||% "none"
  absorbed <- character()
  if (panel %in% c("fixed", "twoway")) absorbed <- c(absorbed, unlist(spec[["unit_constant"]]))
  if (panel == "twoway") absorbed <- c(absorbed, unlist(spec[["time_constant"]]))
  absorbed <- intersect(absorbed, names(g0))
  g <- causal_set_latent(g0, observed, absorbed, absorb = TRUE)
  g_lat <- causal_set_latent(g0, observed, absorbed, absorb = FALSE)
  eff <- spec[["effect"]] %||% "total"
  tot <- causal_choose_set(g, x, y, spec[["adjust"]], "total")
  dir <- if (eff %in% c("direct", "both")) causal_choose_set(g, x, y, spec[["adjust_direct"]] %||% "minimal", "direct") else NULL
  list(g = g, g_latent = g_lat, x = x, y = y, absorbed = absorbed, panel = panel,
       latent = dagitty::latents(g_lat), roles = causal_roles(g, x, y),
       total = tot, direct = dir, effect = eff,
       observed = intersect(names(g0), observed))
}

# ---- independent d-separation (Bayes-ball) for the cross-check --------------
# Reachability-based d-separation (Koller & Friedman 2009, Alg. 3.1); used
# only in dev/causal to verify dagitty's adjustment sets by brute force.
dsep_reachable <- function(pa, x, z) {
  nodes <- names(pa)
  ch <- lapply(stats::setNames(nodes, nodes), function(v) nodes[vapply(nodes, function(w) v %in% pa[[w]], logical(1))])
  anc_z <- z; frontier <- z
  while (length(frontier)) { nf <- setdiff(unlist(pa[frontier]), anc_z); anc_z <- c(anc_z, nf); frontier <- nf }
  visited <- character(); reach <- character()
  todo <- list(list(v = x, dir = "up"))
  while (length(todo)) {
    cur <- todo[[1]]; todo <- todo[-1]
    key <- paste(cur$v, cur$dir)
    if (key %in% visited) next
    visited <- c(visited, key)
    if (!cur$v %in% z) reach <- union(reach, cur$v)
    if (cur$dir == "up" && !cur$v %in% z) {
      for (p in pa[[cur$v]]) todo[[length(todo) + 1]] <- list(v = p, dir = "up")
      for (cc in ch[[cur$v]]) todo[[length(todo) + 1]] <- list(v = cc, dir = "down")
    } else if (cur$dir == "down") {
      if (!cur$v %in% z) for (cc in ch[[cur$v]]) todo[[length(todo) + 1]] <- list(v = cc, dir = "down")
      if (cur$v %in% anc_z) for (p in pa[[cur$v]]) todo[[length(todo) + 1]] <- list(v = p, dir = "up")
    }
  }
  reach
}

# Generalised back-door check of Z for the total effect x -> y (Shpitser,
# VanderWeele & Robins 2010 adjustment criterion), from scratch.
adjustment_criterion <- function(pa, x, y, z) {
  nodes <- names(pa)
  ch <- lapply(stats::setNames(nodes, nodes), function(v) nodes[vapply(nodes, function(w) v %in% pa[[w]], logical(1))])
  desc <- function(v) { d <- v; f <- v; while (length(f)) { nf <- setdiff(unlist(ch[f]), d); d <- c(d, nf); f <- nf }; d }
  anc <- function(v) { d <- v; f <- v; while (length(f)) { nf <- setdiff(unlist(pa[f]), d); d <- c(d, nf); f <- nf }; d }
  cn <- setdiff(intersect(desc(x), anc(y)), x)
  forb <- unique(unlist(lapply(cn, desc)))
  if (length(intersect(z, c(forb, x)))) return(FALSE)
  # proper back-door graph: remove the first edge of every proper causal path
  pa2 <- pa
  for (v in cn) pa2[[v]] <- setdiff(pa2[[v]], x)
  !(y %in% dsep_reachable(pa2, x, z))
}
