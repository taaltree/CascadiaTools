# Real-data reading, data checks, interpretation, figures and methods text
# for the causal module.

fmt <- function(v, d = 3) formatC(signif(v, d), format = "fg", digits = d, flag = "#") |> trimws() |> sub(pattern = "\\.$", replacement = "")
lab_of <- function(sp, v) sp[["labels"]][[v]] %||% v

# ---- read / check -------------------------------------------------------------------
causal_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc$file)) stop("config data: needs `file:` (CSV with one row per site, or per site x time for panels)")
  f <- if (file.exists(dc$file)) dc$file else file.path(config$.dir, dc$file)
  nf <- utils::count.fields(f, sep = ",", quote = "\"", blank.lines.skip = TRUE)
  if (length(unique(nf)) > 1) {
    bad <- which(nf != nf[1])
    stop(sprintf("%s: the header has %d fields but line(s) %s have a different number (%s). Common causes: a decimal comma (12,5), a comma inside an unquoted note, or a missing value with no placeholder. Fix those lines in a copy of the file.",
                 basename(f), nf[1], paste(head(bad + 1L, 10), collapse = ", "), paste(head(nf[bad], 10), collapse = ", ")))
  }
  df <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(dc$columns)) for (node in names(dc$columns)) {
    col <- dc$columns[[node]]
    if (!col %in% names(df)) stop("columns: '", col, "' (for DAG node ", node, ") is not a column of ", basename(f))
    names(df)[names(df) == col] <- node
  }
  spec <- dc[setdiff(names(dc), c("file", "columns"))]
  spec[["dag"]] <- causal_dag_string(dc$dag, dc$edges); spec[["edges"]] <- NULL
  for (k in c("unit_constant", "time_constant", "adjust", "adjust_direct")) if (!is.null(spec[[k]])) spec[[k]] <- unlist(spec[[k]])
  if (is.null(spec[["panel"]])) spec[["panel"]] <- "none"
  list(df = df, spec = spec)
}

causal_check_data <- function(data) {
  pr <- character()
  sp <- data$spec; df <- data$df
  for (k in c("exposure", "outcome", "family")) if (is.null(sp[[k]])) pr <- c(pr, sprintf("config data: `%s:` is required.", k))
  if (length(pr)) return(pr)
  g <- tryCatch(causal_parse(sp[["dag"]]), error = function(e) e)
  if (inherits(g, "error")) return(conditionMessage(g))
  nodes <- names(g)
  for (k in c("exposure", "outcome")) {
    v <- sp[[k]]
    if (!v %in% nodes) pr <- c(pr, sprintf("%s '%s' is not a node of the DAG (nodes: %s).", k, v, paste(nodes, collapse = ", ")))
    if (!v %in% names(df)) pr <- c(pr, sprintf("%s '%s' is not a column of the data (columns: %s). Node names must match column names, or map them with `columns:`.", k, v, paste(names(df), collapse = ", ")))
  }
  if (length(pr)) return(pr)
  if (!sp[["exposure"]] %in% dagitty::ancestors(g, sp[["outcome"]]))
    pr <- c(pr, sprintf("The DAG has no directed path from %s to %s, so it already asserts the causal effect is zero. Add the arrow(s) you want to test.", sp[["exposure"]], sp[["outcome"]]))
  dq <- intersect(dagitty::latents(g), names(df))
  if (length(dq)) pr <- c(pr, sprintf("Node(s) %s are marked [latent] in the DAG but are columns of the data; remove the [latent] tag if they were measured.", paste(dq, collapse = ", ")))
  y <- df[[sp[["outcome"]]]]
  fam <- tryCatch(causal_family_info(sp[["family"]], sp[["link"]]), error = function(e) e)
  if (inherits(fam, "error")) return(c(pr, conditionMessage(fam)))
  yy <- y[!is.na(y)]
  if (!is.numeric(yy)) pr <- c(pr, sprintf("outcome '%s' must be numeric.", sp[["outcome"]]))
  else {
    if (fam$family %in% c("poisson", "negbin") && any(yy < 0 | yy != round(yy))) pr <- c(pr, sprintf("family %s needs non-negative whole-number counts in '%s'.", fam$family, sp[["outcome"]]))
    if (fam$family == "binomial") {
      tr <- if (!is.null(sp[["trials"]])) df[[sp[["trials"]]]] else 1
      if (!is.null(sp[["trials"]]) && is.null(df[[sp[["trials"]]]])) pr <- c(pr, sprintf("trials column '%s' not found.", sp[["trials"]]))
      else if (any(yy < 0 | yy != round(yy)) || any(y > tr, na.rm = TRUE)) pr <- c(pr, sprintf("binomial outcome '%s' must be 0/1 (or successes between 0 and the trials column).", sp[["outcome"]]))
    }
    if (fam$family == "Gamma" && any(yy <= 0)) pr <- c(pr, sprintf("family Gamma needs a strictly positive outcome; '%s' has zeros or negatives.", sp[["outcome"]]))
  }
  xv <- df[[sp[["exposure"]]]]
  if (length(unique(stats::na.omit(xv))) < 2) pr <- c(pr, sprintf("exposure '%s' does not vary.", sp[["exposure"]]))
  obs <- intersect(nodes, names(df))
  unm <- setdiff(nodes, c(names(df), dagitty::latents(g)))
  for (v in unm) {
    near <- names(df)[utils::adist(tolower(v), tolower(names(df))) <= max(2, nchar(v) %/% 5)]
    if (length(near)) pr <- c(pr, sprintf("DAG node '%s' has no column of that name, so it would be treated as UNMEASURED; column(s) %s look similar. Rename the column, map it with `columns:`, or tag the node [latent] if it really was not measured.", v, paste0("'", near, "'", collapse = ", ")))
  }
  for (v in obs) {
    if (is.character(df[[v]])) {
      num <- suppressWarnings(as.numeric(df[[v]]))
      badv <- unique(df[[v]][is.na(num) & !is.na(df[[v]]) & nzchar(trimws(df[[v]]))])
      if (mean(!is.na(num)) > 0.8 && length(badv))
        pr <- c(pr, sprintf("'%s' is numeric except for value(s) %s; fix or blank them (blank = missing).", v, paste0("'", head(badv, 5), "'", collapse = ", ")))
    }
    if (mean(is.na(df[[v]])) > 0.2) pr <- c(pr, sprintf("'%s' is missing in %.0f%% of rows; complete-case analysis would discard too much. Impute or drop it from the DAG with care.", v, 100 * mean(is.na(df[[v]]))))
    if (is.character(df[[v]]) && !v %in% c(sp[["exposure"]]) && length(unique(df[[v]])) > 20 && mean(!is.na(suppressWarnings(as.numeric(df[[v]])))) <= 0.8) pr <- c(pr, sprintf("'%s' is text with %d distinct values; is it really a categorical variable?", v, length(unique(df[[v]]))))
  }
  panel <- sp[["panel"]] %||% "none"
  if (panel != "none") {
    if (is.null(sp[["unit"]]) || !sp[["unit"]] %in% names(df)) pr <- c(pr, "panel: fixed/twoway needs `unit:` naming the site/individual column.")
    else if (all(table(df[[sp[["unit"]]]]) == 1)) pr <- c(pr, sprintf("every %s appears only once: fixed effects need repeated measures of the same units.", sp[["unit"]]))
    else {
      wv <- tapply(df[[sp[["exposure"]]]], df[[sp[["unit"]]]], function(z) length(unique(z)) > 1)
      if (!any(wv)) pr <- c(pr, sprintf("%s never changes within a %s, so a within-%s (fixed-effects) estimate is impossible.", sp[["exposure"]], sp[["unit"]], sp[["unit"]]))
      # misspelt unit ids: ids that differ only in case/spaces, or rare ids
      # one or two characters away from a common id
      u <- as.character(df[[sp[["unit"]]]]); tu <- table(u)
      key <- tolower(gsub("\\s+", "", names(tu)))
      clash <- names(tu)[key %in% key[duplicated(key)]]
      if (length(clash)) pr <- c(pr, sprintf("%s ids differ only in capitals or spaces: %s. They would be treated as different %ss; make them identical.",
                                             sp[["unit"]], paste0("'", head(clash, 6), "'", collapse = ", "), sp[["unit"]]))
      typ <- stats::median(tu)
      rare <- names(tu)[tu <= max(1, typ / 3) & !names(tu) %in% clash]
      common <- setdiff(names(tu), rare)
      for (r in head(rare, 5)) {
        dd <- as.vector(utils::adist(r, common)); near <- if (length(common) && min(dd) <= 2) common[dd == min(dd)] else character()
        if (length(near)) pr <- c(pr, sprintf("%s '%s' has only %d row(s) while most %ss have %d; it looks like a misspelling of %s. Fix the id, or ignore this if it really is a separate %s.",
                                              sp[["unit"]], r, tu[[r]], sp[["unit"]], as.integer(typ), paste0("'", head(near, 3), "'", collapse = ", "), sp[["unit"]]))
      }
      if (!is.null(sp[["time"]]) && sp[["time"]] %in% names(df)) {
        dup <- duplicated(df[, c(sp[["unit"]], sp[["time"]])])
        if (any(dup)) pr <- c(pr, sprintf("%d row(s) repeat a %s x %s combination already in the data (first: %s %s, %s %s). Remove duplicates or check the ids.",
                                          sum(dup), sp[["unit"]], sp[["time"]], sp[["unit"]], u[which(dup)[1]], sp[["time"]], df[[sp[["time"]]]][which(dup)[1]]))
      }
    }
  }
  if (!length(pr)) {
    id <- tryCatch(causal_identify(sp, names(df)), error = function(e) e)
    if (inherits(id, "error")) pr <- c(pr, conditionMessage(id))
    else if (is.null(id$total$set)) pr <- c(pr, causal_not_identified_msg(id))
    else if (identical(id$total$how, "user") && !isTRUE(id$total$valid))
      pr <- c(pr, sprintf("The adjustment set in config {%s} does not block every back-door path (or includes a mediator/collider) in this DAG. Valid minimal sets: %s. Remove `adjust:` to let the module choose.",
                          paste(sp[["adjust"]], collapse = ", "), paste(vapply(id$total$sets, function(s) paste0("{", paste(s, collapse = ", "), "}"), ""), collapse = " ")))
  }
  pr
}

# ---- interpretation -------------------------------------------------------------------
causal_interpret <- function(fit, data) {
  sp <- fit$prep$sp; fam <- fit$prep$fam; x <- fit$prep$x; y <- fit$prep$y
  est <- fit$estimates; nm <- fit$names
  row <- function(p) est[est$param == p, , drop = FALSE]
  xl <- lab_of(sp, x); yl <- sp[["response_label"]] %||% lab_of(sp, y)
  xu <- sp[["unit_text"]][[x]] %||% ""
  a <- row(nm[["ate"]]); pt <- fit$point
  out <- character()
  bin <- fit$prep$ex_type == "binary"
  lv <- sp[["exposure_levels"]] %||% c("0", "1")
  scale_txt <- switch(fam$family, binomial = "probability of", poisson = , negbin = "expected number of", "mean")
  scale_short <- switch(fam$family, binomial = "probability", poisson = , negbin = "count", "response")
  if (bin) {
    s <- sprintf("Average causal effect: if every site had %s = %s instead of %s = %s, the %s %s would change by %s (95%% CI %s to %s): from %s to %s on average.",
                 xl, lv[2], xl, lv[1], scale_txt, yl, fmt(a$estimate), fmt(a$lower), fmt(a$upper), fmt(pt$m0), fmt(pt$m1))
  } else {
    s <- sprintf("Average causal effect: increasing %s by %s%s at every site would change the %s %s by %s (95%% CI %s to %s), from %s to %s on average.",
                 xl, format(fit$prep$contrast), if (nzchar(xu)) paste0(" ", xu) else "", scale_txt, yl,
                 fmt(a$estimate), fmt(a$lower), fmt(a$upper), fmt(pt$m0), fmt(pt$m1))
  }
  if (!is.na(fam$ratio)) {
    r <- row(nm[["ratio"]])
    s <- paste(s, sprintf("As a ratio (%s), %s (%s to %s): a %s%% %s.", gsub("_", " ", fam$ratio), fmt(r$estimate), fmt(r$lower), fmt(r$upper),
                          fmt(abs(100 * (r$estimate - 1)), 2), if (r$estimate >= 1) "increase" else "decrease"))
    s <- sub(" a (8[0-9]?|11|18)(\\.[0-9]+)?% ", " an \\1\\2% ", s)
  }
  if (a$lower > 0 || a$upper < 0) s <- paste(s, "The interval excludes zero, so the direction of the effect is clear under the assumptions below.")
  else s <- paste(s, "The interval includes zero: the data cannot tell whether the effect is positive, negative or absent.")
  out <- c(out, s)
  b <- row(nm[["beta"]])
  per1 <- if (bin) paste("switching", xl, "to", lv[2]) else sprintf("a 1-unit increase in %s%s", xl, if (nzchar(xu)) paste0(" (units: ", xu, ")") else "")
  if (fam$link == "log") out <- c(out, sprintf("On the model scale: %s multiplies the expected %s by %s (95%% CI %s to %s) at every site, because the effect is multiplicative.",
                                              per1, yl, fmt(exp(b$estimate)), fmt(exp(b$lower)), fmt(exp(b$upper))))
  if (fam$link == "logit") out <- c(out, sprintf("On the model scale the conditional odds ratio per %s %s is %s (95%% CI %s to %s). Odds ratios change with which covariates are in the model (non-collapsibility), so use the probability difference above for the causal statement.",
                                                if (bin) "switch to" else paste("1", if (nzchar(xu)) xu else "unit of"), if (bin) lv[2] else xl, fmt(exp(b$estimate)), fmt(exp(b$lower)), fmt(exp(b$upper))))
  if (fam$link == "identity") out <- c(out, sprintf("On the model scale: %s changes %s by %s (95%% CI %s to %s).",
                                                   per1, yl, fmt(b$estimate), fmt(b$lower), fmt(b$upper)))
  d <- row(paste0("direct_ate_", cn_clean(x)))
  if (nrow(d)) {
    med <- names(fit$id$roles)[grepl("^mediator", fit$id$roles)]
    db <- row(paste0("direct_beta_", cn_clean(x)))
    s <- sprintf("Controlled direct effect (holding %s fixed): %s (95%% CI %s to %s) on the %s scale.",
                 paste(vapply(med, function(v) lab_of(sp, v), ""), collapse = ", "), fmt(d$estimate), fmt(d$lower), fmt(d$upper), scale_short)
    if (fam$link %in% c("identity", "log") && sign(db$estimate) == sign(b$estimate) && abs(b$estimate) > 0)
      s <- paste(s, sprintf("On the model scale about %s%% of the total effect runs through the mediator(s) (1 - direct/total = %s), a rough split that ignores the uncertainty in both.",
                            fmt(100 * (1 - db$estimate / b$estimate), 2), fmt(1 - db$estimate / b$estimate, 2)))
    out <- c(out, s)
  }
  wi <- est[grepl("^whatif_", est$param), , drop = FALSE]
  if (nrow(wi)) {
    bad <- wi[!wi$valid_by_dag, ]; good <- wi[wi$valid_by_dag, ]
    pct <- function(v) sprintf("%+.0f%%", 100 * (v - a$estimate) / abs(a$estimate))
    s <- sprintf("Why the adjustment set matters (what-if table): the DAG-based estimate is %s.", fmt(a$estimate))
    if (nrow(bad)) s <- paste(s, "Strategies the DAG says are BIASED give:", paste0(paste(sprintf("%s -> %s (%s)", sub("^what if: ", "", bad$role), fmt(bad$estimate), pct(bad$estimate)), collapse = "; "), "."))
    if (nrow(good)) s <- paste(s, "Other strategies the DAG says are valid give:", paste0(paste(sprintf("%s -> %s", sub("^what if: ", "", good$role), fmt(good$estimate)), collapse = "; "), "."),
                               "Valid strategies should roughly agree with each other; large differences between them suggest the DAG is wrong.")
    out <- c(out, s)
  }
  lat <- fit$id$latent
  pnl <- fit$id$panel_used
  asm <- sprintf(paste0("Causal assumptions, in plain language. These estimates are causal effects ONLY IF: (1) the DAG is right: every common cause of %s and %s is drawn, ",
                        "and each missing arrow really means 'no direct effect'%s; (2) every site could realistically have had the contrasted %s (positivity); ",
                        "(3) one site's %s does not change another site's %s, and the contrast describes a well-defined change (no interference, consistency); ",
                        "(4) the %s model has the right form: %s; (5) the adjusted variables are measured accurately (error in a confounder leaves some confounding behind)%s."),
                 xl, yl, if (length(lat)) paste0("; the unmeasured variables (", paste(lat, collapse = ", "), ") act only as drawn") else "",
                 xl, xl, yl, fam$family,
                 switch(fam$link, log = "the effect is multiplicative and the same in all conditions (no interaction with the adjusted variables)",
                        identity = "the effect is additive and the same in all conditions (no interaction with the adjusted variables)",
                        logit = "effects are linear on the log-odds scale with no interactions"),
                 if (pnl %in% c("fixed", "twoway")) sprintf("; (6) for the fixed-effects estimator: unmeasured confounders do not change over time within a %s%s, and past %s does not feed back into later %s",
                                                            sp[["unit"]], if (pnl == "twoway") " (time-varying ones must be shared by all units each period)" else "", yl, xl) else "")
  out <- c(out, asm)
  if (fam$family %in% c("poisson", "negbin", "binomial"))
    out <- c(out, sprintf(paste0("Detection versus ecology: %s is analysed as recorded, as if detection (capture or sighting efficiency) were the same at every site. ",
                                 "If detectability itself changes with %s, this estimate mixes a change in the ecology with a change in detection; ",
                                 "repeat surveys analysed with an occupancy or N-mixture model are needed to separate the two."), yl, xl))
  ev <- fit$diagnostics$evalue
  if (!is.null(ev)) out <- c(out, sprintf("Sensitivity to unmeasured confounding (E-value): to fully explain away this effect, an unmeasured common cause would have to raise both %s and %s by a factor of %s (risk-ratio scale) beyond the adjusted variables; to make the 95%% CI include no effect, %s. %s",
                                         xl, yl, fmt(as.numeric(ev[["point"]])), if (as.numeric(ev[["ci"]]) <= 1) "no confounding at all would be needed (the CI already includes it)" else paste("a factor of", fmt(as.numeric(ev[["ci"]]))),
                                         "Compare these numbers with the strongest measured confounder you know of."))
  lt <- fit$diagnostics$local_tests
  if (!is.null(lt)) {
    fc <- attr(lt, "fisher_c"); bad <- lt$test[lt$p_holm < 0.05]
    out <- c(out, sprintf("Checking the DAG against the data: the DAG implies %d conditional independencies that can be tested. %s (Fisher's C = %s, df = %d, p = %s). Passing these tests does not prove the DAG: an unmeasured confounder of %s and %s leaves no trace in them.",
                          nrow(lt), if (length(bad)) paste0("The data contradict ", length(bad), " of them (", paste(bad, collapse = "; "), "): the DAG is missing arrows or common causes, so revise it") else "None is contradicted by the data",
                          fmt(fc[["C"]]), as.integer(fc[["df"]]), fmt(fc[["p"]], 2), xl, yl))
  }
  dp <- fit$diagnostics$dharma_p
  if (!is.null(dp)) out <- c(out, if (any(dp < 0.05)) sprintf("Model check: the %s outcome model shows problems (DHARMa %s p < 0.05); the effect estimate may be off.%s", fam$family, paste(names(dp)[dp < 0.05], collapse = ", "),
                                                    if (fam$family == "poisson" && dp[["dispersion"]] < 0.05) " Counts are more variable than Poisson: use family negbin." else " Check the residual plot; a single p < 0.05 among several checks also happens by chance in 1 of 20 good models.")
                               else sprintf("Model check: simulated residuals of the outcome model show no problems (DHARMa uniformity p = %s, dispersion p = %s).", fmt(dp[["uniformity"]], 2), fmt(dp[["dispersion"]], 2)))
  out <- c(out, paste0("Limits from simulation testing (VALIDATION.md): with the correct DAG the estimates recovered the true effect (bias within 5%, 95% intervals covering the truth 92-97% of the time) in every tested design with 100 or more sites or 20 or more sites x 5 years; ",
                       "adjusting for a mediator or collider, leaving out a confounder, or ignoring site/year confounding in panels gave estimates that missed the truth by 18-320%. ",
                       "With only 40 sites the response-scale effect on counts was biased upward by about 18% (ratios and model-scale coefficients were not)",
                       if (nrow(data$df) < 80 && fam$family %in% c("poisson", "negbin")) sprintf("; this dataset has %d rows, so read the count-scale effect with caution and rely on the ratio.", nrow(data$df)) else ".",
                       " None of this testing can check the DAG itself."))
  out
}

# ---- figures ---------------------------------------------------------------------------
# Publication suite (docs/FIGURES.md): dag, effect, adjustment_comparison,
# implied_independencies, residuals. Every label goes through fig_label():
#   - DAG nodes / variables: their column name (e.g. canopy_pct)
#   - legend entries: role_exposure, role_outcome, role_adjusted, role_absorbed,
#     role_unmeasured, role_not_adjusted, path_causal, path_other,
#     whatif_main, whatif_valid, whatif_biased, test_consistent, test_contradicted
#   - axes: effect_x, effect_y, comparison_x, tests_x, tests_y, resid_x, resid_y
# Theme, base size and journal width come from the core (apply_fig_opts).

causal_cap1 <- function(s) if (!nzchar(s)) s else paste0(toupper(substring(s, 1, 1)), substring(s, 2))
# lower-case the first letter for use mid-sentence, unless it starts an acronym (e.g. "SST")
causal_decap <- function(s) ifelse(grepl("^[A-Z][a-z]", s), paste0(tolower(substring(s, 1, 1)), substring(s, 2)), s)
# label for use inside running text: figure labels (written for axes, so
# capitalised) are decapitalised; data labels are used as written
causal_vlab_text <- function(v, sp, opts) if (!is.null(opts$labels[[v]])) causal_decap(opts$labels[[v]]) else causal_vlab(v, sp, opts)
causal_vlab <- function(v, sp, opts) fig_label(v, if (is.null(sp[["labels"]][[v]])) gsub("_", " ", v) else lab_of(sp, v), opts)

# Replace DAG node names inside free text (strategy names, test strings) with
# their labels. Two passes through placeholders so a label that contains
# another node's name is not rewritten twice.
causal_relabel_text <- function(txt, vars, sp, opts) {
  vars <- unique(vars[nzchar(vars)]); vars <- vars[order(-nchar(vars))]
  esc <- function(s) paste0("\\Q", s, "\\E")
  for (k in seq_along(vars))
    txt <- gsub(paste0("(?<![A-Za-z0-9_.])", esc(vars[k]), "(?![A-Za-z0-9_])"), sprintf("\001%d\002", k), txt, perl = TRUE)
  for (k in seq_along(vars)) txt <- gsub(sprintf("\001%d\002", k), causal_vlab_text(vars[k], sp, opts), txt, fixed = TRUE)
  txt
}

# Legend position. These figures read best with the legend underneath, so the
# core default ("right") becomes "bottom"; any explicit per-figure choice and
# any other global choice (top, none, ...) is respected. Set via the guide
# itself because the core replaces the module's theme.
causal_legend_pos <- function(opts, id) {
  pf <- opts$per_figure[[id]]$legend
  if (!is.null(pf)) return(pf)
  lp <- opts$legend %||% "right"
  if (identical(lp, "right")) "bottom" else lp
}
causal_guide <- function(lp, ...) if (identical(lp, "none")) "none" else ggplot2::guide_legend(position = lp, ...)

causal_ci_lab <- function(level) sprintf("%s%%", trimws(formatC(100 * level, format = "fg", digits = 3)))

causal_layout <- function(g, edges) {
  co <- dagitty::coordinates(g)
  if (length(co$x) && !anyNA(co$x)) return(data.frame(name = names(co$x), x = co$x, y = -co$y))
  v <- names(g); lvl <- stats::setNames(rep(0, length(v)), v)
  for (i in seq_along(v)) for (k in seq_len(nrow(edges))) lvl[edges$to[k]] <- max(lvl[edges$to[k]], lvl[edges$from[k]] + 1)
  # Sugiyama layering (igraph) orders nodes within layers to reduce crossings
  ly <- tryCatch({
    ig <- igraph::graph_from_data_frame(edges, vertices = data.frame(name = v))
    L <- igraph::layout_with_sugiyama(ig, layers = unname(lvl[igraph::V(ig)$name]))$layout
    data.frame(name = igraph::V(ig)$name, x = L[, 2], y = L[, 1])
  }, error = function(e) NULL)
  if (!is.null(ly)) {
    # layers run left -> right; spread nodes within a layer
    ly$x <- max(ly$x) - ly$x
    ly$y <- stats::ave(ly$y, ly$x, FUN = function(z) rank(z) - (length(z) + 1) / 2)
    return(ly)
  }
  nodes <- data.frame(name = v, x = lvl[v])
  nodes$y <- stats::ave(nodes$x, nodes$x, FUN = function(z) seq_along(z) - (length(z) + 1) / 2)
  nodes
}

causal_fig_dag <- function(fit, opts) {
  g <- causal_parse(fit$dag_string)
  e <- causal_edges(g)
  nodes <- causal_layout(g, e)
  nodes$x <- nodes$x * 1.6
  sp <- fit$prep$sp
  x <- fit$prep$x; y <- fit$prep$y
  adj <- fit$sets$chosen
  rl <- c(exposure = fig_label("role_exposure", "Exposure", opts), outcome = fig_label("role_outcome", "Outcome", opts),
          adjusted = fig_label("role_adjusted", "Adjusted for", opts), absorbed = fig_label("role_absorbed", "Absorbed by fixed effects", opts),
          unmeasured = fig_label("role_unmeasured", "Unmeasured", opts), not_adjusted = fig_label("role_not_adjusted", "Not adjusted", opts))
  role <- ifelse(nodes$name == x, "exposure", ifelse(nodes$name == y, "outcome",
          ifelse(nodes$name %in% adj & nodes$name %in% fit$id$absorbed, "absorbed",
          ifelse(nodes$name %in% adj, "adjusted", ifelse(nodes$name %in% fit$id$latent, "unmeasured", "not_adjusted")))))
  nodes$role <- factor(rl[role], levels = unname(rl))
  # second encoding for greyscale: exposure/outcome bold, unmeasured italic
  nodes$face <- ifelse(role %in% c("exposure", "outcome"), "bold", ifelse(role == "unmeasured", "italic", "plain"))
  nodes$label <- vapply(nodes$name, function(v) paste(strwrap(causal_vlab(v, sp, opts), 14), collapse = "\n"), "")
  pl <- c(causal = fig_label("path_causal", "Causal path", opts), other = fig_label("path_other", "Other arrow", opts))
  cn <- causal_cn(g, x, y)
  e$path <- ifelse(e$from %in% c(x, cn) & e$to %in% cn, pl[["causal"]], pl[["other"]])
  e <- merge(merge(e, stats::setNames(nodes[, c("name", "x", "y")], c("from", "x0", "y0"))),
             stats::setNames(nodes[, c("name", "x", "y")], c("to", "x1", "y1")))
  # Lines run centre to centre underneath the (opaque) labels; the arrowhead
  # sits part-way along each edge so it is never hidden. An edge that would
  # pass through another node is bent around it.
  segs <- list()
  for (k in seq_len(nrow(e))) {
    p0 <- c(e$x0[k], e$y0[k]); p1 <- c(e$x1[k], e$y1[k]); d <- p1 - p0; L <- sqrt(sum(d^2))
    others <- nodes[!nodes$name %in% c(e$from[k], e$to[k]), c("x", "y")]
    tt <- ((others$x - p0[1]) * d[1] + (others$y - p0[2]) * d[2]) / L^2
    dist <- abs((others$x - p0[1]) * d[2] - (others$y - p0[2]) * d[1]) / L
    blocked <- any(tt > 0.05 & tt < 0.95 & dist < 0.32)
    if (blocked) {
      perp <- c(-d[2], d[1]) / L
      bend <- (p0 + p1) / 2 + perp * 0.65 * (if (k %% 2) 1 else -1)
      segs[[length(segs) + 1]] <- data.frame(x = p0[1], y = p0[2], xend = bend[1], yend = bend[2], arrow = TRUE, path = e$path[k])
      segs[[length(segs) + 1]] <- data.frame(x = bend[1], y = bend[2], xend = p1[1], yend = p1[2], arrow = FALSE, path = e$path[k])
    } else {
      mid <- p0 + 0.58 * d
      segs[[length(segs) + 1]] <- data.frame(x = p0[1], y = p0[2], xend = mid[1], yend = mid[2], arrow = TRUE, path = e$path[k])
      segs[[length(segs) + 1]] <- data.frame(x = mid[1], y = mid[2], xend = p1[1], yend = p1[2], arrow = FALSE, path = e$path[k])
    }
  }
  sg <- do.call(rbind, segs)
  # Arrowheads are drawn as short solid stubs at the end of each arrow
  # segment, so heads stay crisp on dashed edges (a dashed head is a smudge).
  hd <- sg[sg$arrow, ]
  ln <- sqrt((hd$xend - hd$x)^2 + (hd$yend - hd$y)^2)
  hd$x <- hd$xend - 0.04 * (hd$xend - hd$x) / pmax(ln, 1e-9)
  hd$y <- hd$yend - 0.04 * (hd$yend - hd$y) / pmax(ln, 1e-9)
  # Node text is in mm, the layout in data units: shrink the labels when the
  # panel is narrow (single/onehalf width, or a legend on the side) so nodes
  # neither overlap nor get clipped at the panel edge.
  lp <- causal_legend_pos(opts, "dag")
  wmm <- fig_width_mm(opts$per_figure$dag$width %||% "double", opts) %||% 170
  pan <- wmm - 8 - (if (lp %in% c("left", "right")) 48 else 0)
  mpu <- pan / (diff(range(nodes$x)) + 1.4)
  # keep text >= ~6 pt (2.1 mm) and wrap names tighter instead when space is short
  avail <- 0.85 * 1.6 * mpu
  lab_size <- max(2.1, min(2.3, 2.3 * avail / 21))
  wrap <- max(8, min(14, floor(avail / (0.62 * lab_size)) - 1))
  nodes$label <- vapply(nodes$name, function(v) paste(strwrap(causal_vlab(v, sp, opts), wrap), collapse = "\n"), "")
  # Role colours are fixed per entity (exposure always orange, outcome sky blue,
  # adjustment set bluish green, the same green as the main estimate in
  # adjustment_comparison). Greyscale: graded greys + font face.
  cols <- if (isTRUE(opts$greyscale))
    c(exposure = "#8c8c8c", outcome = "#b3b3b3", adjusted = "#d9d9d9", absorbed = "#ececec", unmeasured = "#f7f7f7", not_adjusted = "#FFFFFF")
  else c(exposure = "#E69F00", outcome = "#56B4E9", adjusted = "#009E73", absorbed = "#CC79A7", unmeasured = "#BBBBBB", not_adjusted = "#FFFFFF")
  names(cols) <- rl[names(cols)]
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = sg, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linetype = path, linewidth = path),
                          colour = eco_pal$ink, key_glyph = "path") +
    ggplot2::geom_segment(data = hd, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = path), linetype = "solid",
                          colour = eco_pal$ink, arrow = ggplot2::arrow(length = ggplot2::unit(1.8, "mm"), type = "closed"),
                          show.legend = FALSE) +
    ggplot2::geom_label(data = nodes, ggplot2::aes(x, y, label = label, fill = role, fontface = face), size = lab_size, linewidth = 0.25,
                        lineheight = 0.9, label.padding = ggplot2::unit(0.15, "lines")) +
    ggplot2::scale_fill_manual(values = cols[levels(droplevels(nodes$role))], name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22"), unname(pl)), name = NULL) +
    ggplot2::scale_linewidth_manual(values = stats::setNames(c(0.7, 0.35), unname(pl)), name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.7)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::guides(x = "none", y = "none",
                    fill = causal_guide(lp, nrow = 2, order = 1, override.aes = list(label = "  ")),
                    linetype = causal_guide(lp, nrow = 2, order = 2), linewidth = causal_guide(lp, nrow = 2, order = 2)) +
    ggplot2::labs(x = NULL, y = NULL)
  attr(p, "fig_size") <- list(width = "double", height = 45 + 20 * max(table(nodes$x)))
  p
}

causal_fig_effect <- function(fit, data, opts) {
  m <- fit$object; pr <- fit$prep; sp <- pr$sp; x <- pr$x
  df <- fit$data_used
  nd <- df; if (!is.null(sp[["offset"]])) nd$.log_offset <- log(sp[["offset_per"]] %||% 1)
  pnl <- fit$id$panel_used
  V <- tryCatch(causal_vcov(m, list(df = df, sp = sp), pnl, fit$settings), error = function(e) stats::vcov(m))
  ci <- opts$ci %||% 0.95
  yl <- fig_label("effect_y", causal_cap1(fig_label(pr$y, sp[["response_label"]] %||% lab_of(sp, pr$y), opts)), opts)
  xu <- sp[["unit_text"]][[x]] %||% ""
  xl0 <- causal_cap1(causal_vlab(x, sp, opts))
  xl <- fig_label("effect_x", if (nzchar(xu) && pr$ex_type != "binary" && is.null(opts$labels[[x]])) sprintf("%s (%s)", xl0, xu) else xl0, opts)
  # long axis titles wrap so they are not clipped at single-column width
  yl <- paste(strwrap(yl, 30), collapse = "\n"); xl <- paste(strwrap(xl, 50), collapse = "\n")
  if (pr$ex_type == "binary") {
    ap <- as.data.frame(marginaleffects::avg_predictions(m, variables = stats::setNames(list(c(0, 1)), x), newdata = nd, vcov = V, by = x, conf_level = ci))
    lv <- sp[["exposure_levels"]] %||% c("0", "1")
    lv <- vapply(as.character(lv), function(l) fig_label(paste0(x, "_", l), l, opts), "")
    pred <- data.frame(level = factor(lv[ap[[x]] + 1], levels = lv), estimate = ap$estimate, lower = ap$conf.low, upper = ap$conf.high)
    p <- fig_cat_effect(pred, xl, yl)
    if (pr$fam$family %in% c("poisson", "negbin", "binomial")) p <- p + ggplot2::expand_limits(y = 0)
    return(p)
  }
  grid <- seq(stats::quantile(df[[x]], 0.02), stats::quantile(df[[x]], 0.98), length.out = 60)
  ap <- as.data.frame(marginaleffects::avg_predictions(m, variables = stats::setNames(list(grid), x), newdata = nd, vcov = V, by = x, conf_level = ci))
  # x in original units: the module fits on the user's own scale, so fig_x is
  # a no-op unless the user supplies `units:` or the data record a scaling.
  pred <- data.frame(x = fig_x(ap[[x]], x, opts, data$scaling), estimate = ap$estimate, lower = ap$conf.low, upper = ap$conf.high)
  obs <- if (isTRUE(opts$show_data)) data.frame(x = fig_x(df[[x]], x, opts, data$scaling)) else NULL
  prob <- pr$fam$family == "binomial" && is.null(sp[["trials"]])
  p <- fig_cond_effect(pred, xl, yl, obs = obs, ylim01 = prob) +
    # extra room on the right so the last tick label is not clipped
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.10)))
  if (!prob) {
    p <- p + ggplot2::scale_y_continuous(breaks = scales::breaks_extended(6))
    # counts / positive responses: show the effect against zero so the
    # relative size is not exaggerated
    if (pr$fam$family %in% c("poisson", "negbin", "binomial") || isTRUE(all(df[[pr$y]] >= 0)))
      p <- p + ggplot2::expand_limits(y = 0)
  }
  p
}

# Wald interval at another level from the stored SD (the fit stores 95%).
causal_rescale_ci <- function(d, level) {
  if (abs(level - 0.95) < 1e-9 || !"sd" %in% names(d)) return(d)
  z <- stats::qnorm(1 - (1 - level) / 2)
  ok <- is.finite(d$sd)
  d$lower[ok] <- d$estimate[ok] - z * d$sd[ok]; d$upper[ok] <- d$estimate[ok] + z * d$sd[ok]
  d
}

causal_fig_comparison <- function(fit, opts) {
  est <- fit$estimates; nm <- fit$names; sp <- fit$prep$sp
  ci <- opts$ci %||% 0.95
  main <- est[est$param == nm[["ate"]], ]
  wi <- est[grepl("^whatif_", est$param), ]
  if (!nrow(wi)) return(NULL)
  sl <- c(main = fig_label("whatif_main", "DAG set (main estimate)", opts), valid = fig_label("whatif_valid", "Valid by the DAG", opts),
          biased = fig_label("whatif_biased", "Biased by the DAG", opts))
  d <- rbind(transform(main, status = sl[["main"]]),
             transform(wi, status = ifelse(wi$valid_by_dag, sl[["valid"]], sl[["biased"]])))
  d <- causal_rescale_ci(d, ci)
  d$status <- factor(d$status, levels = unname(sl))
  vars <- names(causal_parse(fit$dag_string))
  d$label <- vapply(causal_relabel_text(d$strategy, vars, sp, opts), function(z) paste(strwrap(z, 34), collapse = "\n"), "")
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  cols <- if (isTRUE(opts$greyscale)) c("#1a1a1a", "#666666", "#a6a6a6") else c("#009E73", "#0072B2", "#D55E00")
  names(cols) <- sl
  shp <- stats::setNames(c(16, 17, 15), sl)
  lp <- causal_legend_pos(opts, "adjustment_comparison")
  xu <- sp[["unit_text"]][[fit$prep$x]] %||% ""
  what <- if (fit$prep$ex_type == "binary") "switching" else paste0("+", format(fit$prep$contrast), if (nzchar(xu)) paste0(" ", xu, " of") else "")
  xl <- fig_label("comparison_x", sprintf("Average causal effect of %s %s\non %s (%s CI)", what, causal_vlab_text(fit$prep$x, sp, opts),
                                          (if (is.null(opts$labels[[fit$prep$y]])) sp[["response_label"]] %||% lab_of(sp, fit$prep$y) else causal_decap(opts$labels[[fit$prep$y]])), causal_ci_lab(ci)), opts)
  p <- ggplot2::ggplot(d, ggplot2::aes(estimate, label, colour = status, shape = status)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = eco_pal$null) +
    ggplot2::geom_vline(xintercept = main$estimate, colour = cols[[1]], alpha = 0.6) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.3, linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
    ggplot2::scale_shape_manual(values = shp, name = NULL, drop = TRUE) +
    ggplot2::guides(colour = causal_guide(lp), shape = causal_guide(lp)) +
    ggplot2::labs(x = xl, y = NULL)
  attr(p, "fig_size") <- list(width = "double", height = 38 + 7.5 * nrow(d))
  p
}

causal_fig_tests <- function(fit, opts) {
  lt <- fit$diagnostics$local_tests
  if (is.null(lt) || !nrow(lt) || !"lower" %in% names(lt) || all(is.na(lt$lower))) return(NULL)
  sp <- fit$prep$sp
  ci <- opts$ci %||% 0.95
  # dagitty's "cis" intervals are Fisher-z intervals; recompute them at
  # opts$ci when we can reproduce the stored 95% bounds, otherwise keep 95%.
  lvl <- 0.95
  n <- nrow(fit$data_used)
  nz <- vapply(strsplit(gsub("_\\|\\|_", "~", lt$test), " \\| "), function(s) if (length(s) < 2) 0L else length(strsplit(s[2], ",")[[1]]), 1L)
  se <- 1 / sqrt(pmax(n - nz - 3, 1))
  chk <- tanh(atanh(lt$estimate) - stats::qnorm(0.975) * se)
  if (abs(ci - 0.95) > 1e-9 && isTRUE(max(abs(chk - lt$lower), na.rm = TRUE) < 0.005)) {
    z <- stats::qnorm(1 - (1 - ci) / 2)
    lt$lower <- tanh(atanh(lt$estimate) - z * se); lt$upper <- tanh(atanh(lt$estimate) + z * se); lvl <- ci
  }
  vars <- names(causal_parse(fit$dag_string))
  lab <- gsub(" _\\|\\|_ ", " vs ", causal_relabel_text(lt$test, vars, sp, opts))
  lt$label <- vapply(lab, function(z) paste(strwrap(z, 62, exdent = 4), collapse = "\n"), "")
  lt$label <- factor(lt$label, levels = rev(lt$label))
  sl <- c(ok = fig_label("test_consistent", "Consistent with the DAG", opts),
          bad = fig_label("test_contradicted", "Contradicted (Holm p < 0.05)", opts))
  lt$status <- factor(ifelse(lt$p_holm < 0.05, sl[["bad"]], sl[["ok"]]), levels = unname(sl))
  cols <- stats::setNames(if (isTRUE(opts$greyscale)) c("#666666", "#1a1a1a") else c("#0072B2", "#D55E00"), sl)
  lp <- causal_legend_pos(opts, "implied_independencies")
  p <- ggplot2::ggplot(lt, ggplot2::aes(estimate, label, colour = status, shape = status)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = eco_pal$null) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.3, linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 15), sl), name = NULL, drop = TRUE) +
    ggplot2::guides(colour = causal_guide(lp), shape = causal_guide(lp)) +
    ggplot2::labs(x = fig_label("tests_x", sprintf("Partial correlation the DAG says is 0 (%s CI)", causal_ci_lab(lvl)), opts),
                  y = fig_label("tests_y", "Independence tested (A vs B | given)", opts))
  nl <- sum(lengths(regmatches(as.character(lt$label), gregexpr("\n", as.character(lt$label))))) + nrow(lt)
  attr(p, "fig_size") <- list(width = "double", height = 34 + 5.6 * nl)
  attr(p, "ci_level") <- lvl
  p
}

causal_fig_residuals <- function(fit, opts = fig_opts()) {
  r <- fit$dharma$res
  if (is.null(r)) return(NULL)
  s <- sort(r$scaledResiduals); n <- length(s)
  d <- data.frame(expected = (seq_len(n) - 0.5) / n, observed = s)
  p <- ggplot2::ggplot(d, ggplot2::aes(expected, observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = eco_pal$null, linetype = 2) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6, colour = eco_pal$ink) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = fig_label("resid_x", "Expected quantile (uniform)", opts),
                  y = fig_label("resid_y", "Simulated (DHARMa) residual", opts))
  attr(p, "fig_size") <- list(aspect = 0.85)
  p
}

causal_figures <- function(fit, data, opts = fig_opts()) {
  pr <- fit$prep; sp <- pr$sp
  ci <- opts$ci %||% 0.95
  cil <- causal_ci_lab(ci)
  xl <- causal_vlab_text(pr$x, sp, opts)
  yl <- if (is.null(opts$labels[[pr$y]])) sp[["response_label"]] %||% lab_of(sp, pr$y) else causal_decap(opts$labels[[pr$y]])
  safe <- function(expr) tryCatch(expr, error = function(e) { message("causal figure skipped: ", conditionMessage(e)); NULL })
  cap <- function(p, txt) if (is.null(p)) NULL else with_caption(p, txt)
  adj <- fit$sets$chosen
  adj_txt <- if (length(adj)) paste(vapply(adj, causal_vlab_text, "", sp = sp, opts = opts), collapse = ", ") else "nothing"
  tests <- safe(causal_fig_tests(fit, opts))
  tl <- if (is.null(tests)) "95%" else causal_ci_lab(attr(tests, "ci_level"))
  gs <- isTRUE(opts$greyscale)
  cw <- function(colour, grey) if (gs) grey else colour
  out <- list(
    dag = cap(safe(causal_fig_dag(fit, opts)),
      paste0(sprintf("Causal diagram (DAG) assumed for the analysis. Exposure: %s (%s); outcome: %s (%s). %s nodes are the adjustment set used for the total effect (%s)",
              xl, cw("orange", "dark grey, bold"), yl, cw("blue", "mid grey, bold"), cw("Green", "Light grey"), adj_txt),
             # mention only the roles actually drawn
             if (any(adj %in% fit$id$absorbed)) sprintf("; %s nodes are absorbed by fixed effects", cw("pink", "very light grey")) else "",
             if (length(setdiff(fit$id$latent, adj))) sprintf("; %s nodes are unmeasured", cw("grey", "italic")) else "",
             if (length(setdiff(names(causal_parse(fit$dag_string)), c(pr$x, pr$y, adj, fit$id$latent)))) "; white nodes are not adjusted for" else "",
             ". Solid arrows lie on causal paths from the exposure to the outcome; dashed arrows are other assumed direct effects. Each missing arrow is an assumption of no direct effect.")),
    effect = cap(safe(causal_fig_effect(fit, data, opts)),
      if (pr$ex_type == "binary") sprintf("Average %s if every site had each level of %s (g-computation over the observed sites, adjusted for %s), with %s confidence intervals.", yl, xl, adj_txt, cil)
      else sprintf("Causal dose-response curve: average predicted %s if every site were set to each value of %s, keeping each site's own values of the adjustment variables (%s; g-computation). Shaded band: %s confidence interval (conditional on the sites sampled). Rug: observed values of %s. The curve spans the 2nd-98th percentiles of the observed exposure.",
                   yl, xl, adj_txt, cil, xl)),
    adjustment_comparison = cap(safe(causal_fig_comparison(fit, opts)),
      sprintf("What-if comparison: the same average causal effect estimated under different adjustment choices, with %s confidence intervals. %s circle and vertical line: the DAG's adjustment set (main estimate). %s triangles: other sets the DAG says are valid. %s squares: choices the DAG says are biased (omitting a confounder, adjusting for a mediator or a collider). Dashed line: no effect.", cil, cw("Green", "Black"), cw("Blue", "Dark grey"), cw("Orange", "Light grey"))),
    implied_independencies = cap(tests,
      sprintf("Testable implications of the DAG: partial correlations that the DAG says should be zero (A vs B given the variables after |), with %s confidence intervals (dagitty::localTests). %s Consistency does not prove the DAG is right.", tl,
              if (is.null(tests)) "" else if (any(fit$diagnostics$local_tests$p_holm < 0.05, na.rm = TRUE))
                sprintf("%s squares are contradicted by the data after Holm correction; %s circles are consistent with the DAG.", cw("Orange", "Black"), cw("blue", "grey"))
              else sprintf("All %d are consistent with the DAG (none contradicted after Holm correction).", nrow(fit$diagnostics$local_tests)))),
    residuals = cap(safe(causal_fig_residuals(fit, opts)),
      "Outcome model check: DHARMa simulated quantile residuals of the outcome model plotted against uniform quantiles. Points should follow the dashed 1:1 line; S-shapes indicate over- or underdispersion.")
  )
  Filter(function(p) !is.null(p) && inherits(p, "ggplot"), out)
}

causal_plots <- function(fit, data) {
  f <- causal_figures(fit, data, fig_opts())
  lapply(f, function(p) p + theme_eco(11))
}

# ---- methods text -------------------------------------------------------------------
causal_methods_text <- function(fit, data, settings) {
  pr <- fit$prep; sp <- pr$sp; fam <- pr$fam; s <- fit$sets
  st <- merge_list(CAUSAL_DEFAULT_SETTINGS, settings)
  pnl <- fit$id$panel_used
  paste0(
    sprintf("We encoded our assumptions about the causal structure as a directed acyclic graph (DAG; %d variables, Pearl 2009; Arif & MacNeil 2023) and used the back-door/adjustment criterion (dagitty %s; Textor et al. 2016) to find the minimal sets of measured covariates sufficient to identify the total causal effect of %s on %s. ",
            length(names(causal_parse(fit$dag_string))), as.character(utils::packageVersion("dagitty")), pr$x, pr$y),
    sprintf("There %s %d such set%s; we adjusted for {%s}%s. ", if (length(s$minimal) == 1) "was" else "were", length(s$minimal), if (length(s$minimal) == 1) "" else "s",
            paste(setdiff(s$chosen, fit$id$absorbed), collapse = ", "), if (s$how == "optimal") " (the optimal adjustment set; Henckel et al. 2022)" else if (s$how == "user") " (set specified a priori)" else ""),
    if (pnl %in% c("fixed", "twoway")) sprintf("Because %s was unmeasured but assumed constant within each %s, we used a within-%s (fixed-effects) estimator with %s indicators%s and cluster-robust standard errors by %s (Byrnes & Dee 2025). ",
                                               paste(fit$id$absorbed, collapse = ", "), sp[["unit"]], sp[["unit"]], sp[["unit"]], if (pnl == "twoway") paste0(" and ", sp[["time"]], " indicators (two-way fixed effects)") else "", sp[["unit"]]) else "",
    sprintf("The outcome was modelled with a %s GLM (%s link) containing the exposure and the adjustment set (R %s, %s). ", fam$family, fam$link, paste(R.version$major, R.version$minor, sep = "."),
            if (fam$family == "negbin") "MASS::glm.nb" else "stats::glm"),
    sprintf("The average causal effect of %s on the response scale was obtained by g-computation (standardisation): predictions for every site with its exposure %s were averaged and contrasted (marginaleffects %s; Arel-Bundock et al. 2024). ",
            pr$x, if (pr$ex_type == "binary") "set to each level" else paste("increased by", format(pr$contrast), "versus as observed"),
            as.character(utils::packageVersion("marginaleffects"))),
    switch(st$ci_method,
           unconditional = "Standard errors combine the delta-method variance of the model coefficients with the sampling variance of the covariate distribution (influence-function term), so intervals refer to the population of sites rather than only the sites sampled. ",
           delta = "Standard errors use the delta method, conditional on the sites sampled. ",
           bootstrap = sprintf("Confidence intervals are percentile intervals from %d nonparametric bootstrap resamples%s. ", st$n_boot, if (pnl != "none") paste0(" of whole ", sp[["unit"]], "s") else "")),
    if (fit$id$effect %in% c("direct", "both")) "The controlled direct effect was estimated in a second model with the covariate set dagitty identifies for the direct effect (which includes the mediators). " else "",
    "We tested the DAG's implied conditional independencies against the data (dagitty::localTests; one test per non-adjacent pair, combined with Fisher's C; Shipley 2000), checked the outcome model with simulated residuals (DHARMa), and computed E-values for sensitivity to unmeasured confounding (VanderWeele & Ding 2017). ",
    "For comparison we also estimated the effect under alternative adjustment strategies (no adjustment, omitting each confounder, adding mediators, colliders or all measured variables; Cinelli et al. 2022). ",
    "The procedure was validated by simulation from structural causal models with known effects (ecotools causal module, VALIDATION.md)."
  )
}
