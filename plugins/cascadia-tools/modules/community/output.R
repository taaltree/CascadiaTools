# Interpretation, figures and methods text for the community module.

com_dist_name <- function(d) switch(d, bray = "Bray-Curtis", jaccard = "Jaccard", euclidean = "Euclidean", horn = "Horn-Morisita",
                                  morisita = "Morisita", kulczynski = "Kulczynski", gower = "Gower", d)

com_get <- function(est, p, col = "estimate") { v <- est[[col]][est$param == p]; if (length(v)) v[1] else NA_real_ }
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
fmt_p <- function(p) vapply(p, function(q) if (is.na(q)) "not tested" else if (q < 0.001) "p < 0.001" else sprintf("p = %.3f", q), "")
com_lab <- function(key, opts = NULL, data = NULL) {
  def <- if (!is.null(data$labels[[key]])) data$labels[[key]] else gsub("_", " ", key)
  fig_label(key, def, opts)
}

com_interpret <- function(fit, data) {
  est <- fit$estimates; o <- fit$object; out <- character()
  part <- names(o$partition %||% data$partition)
  lv <- function(k) com_lab(k, NULL, data)
  ci <- function(p) sprintf("%s (95%% CI %s to %s)", pct(com_get(est, p)), pct(com_get(est, p, "lower")), pct(com_get(est, p, "upper")))
  adjtxt <- if (identical(o$settings$r2, "unadjusted")) "unadjusted R2 (biased upward; for teaching only)" else "adjusted R2"
  tr <- o$settings$transform
  if (length(part)) {
    tot <- com_get(est, "R2adj_env_total")
    s <- sprintf("Community composition (%s-transformed %s, redundancy analysis) was associated with the measured environment: together the variables accounted for %s of the among-site variation (%s, %s), leaving %s unexplained.",
                 tr, if (isTRUE(o$presence)) "presence/absence" else "abundances", ci("R2adj_env_total"), adjtxt,
                 fmt_p(com_get(est, "R2adj_env_total", "p_value")), pct(1 - tot))
    out <- c(out, s)
    if (length(part) > 1) {
      seg <- vapply(part, function(k) sprintf("%s alone %s (%s); %s after accounting for the other variables %s (%s)",
                                              lv(k), ci(paste0("R2adj_", k)), fmt_p(com_get(est, paste0("R2adj_", k), "p_value")),
                                              tolower(lv(k)), ci(paste0("R2adj_unique_", k)),
                                              fmt_p(com_get(est, paste0("R2adj_unique_", k), "p_value"))), "")
      out <- c(out, paste0("Variation partitioning: ", paste(seg, collapse = "; "), sprintf(
        ". The shared fraction, %s, is variation that cannot be attributed to one variable because the variables are correlated with each other; it is not an interaction and has no test.", ci("R2adj_shared"))))
      sh <- com_get(est, "R2adj_shared")
      if (!is.na(sh) && sh < 0) out <- c(out, "The shared fraction is negative, which happens when variables act in opposite directions (suppression) or by chance when it is near zero; read it as 'about zero'.")
    }
    big <- part[which.max(vapply(part, function(k) com_get(est, paste0("R2adj_", k)), 0))]
    out <- c(out, sprintf("These percentages are shares of variation in community composition, not effect sizes for any single species: %s explaining %s means sites with similar %s tend to share species and their relative abundances. The association is correlational: unmeasured variables that co-vary with %s could produce the same pattern.",
                          lv(big), pct(com_get(est, paste0("R2adj_", big))), tolower(lv(big)), tolower(lv(big))))
  }
  dl <- com_get(est, "dca_axis1_length_sd")
  if (!is.na(dl)) out <- c(out, sprintf(
    "The first DCA axis is %.1f SD units long. %s", dl,
    if (dl > 4) "Above about 4 SD most species are replaced along the main gradient (unimodal responses). RDA on raw abundances would then be inappropriate; the Hellinger (or chord) transformation used here makes RDA suitable for such data (Legendre & Gallagher 2001), and CCA (unimodal, chi-square distance) is reported below as a check."
    else if (dl < 3) "Below about 3 SD species respond roughly linearly over the sampled range, so RDA is appropriate."
    else "Between 3 and 4 SD both linear (RDA) and unimodal (CCA) methods are reasonable; compare them."))
  cc <- com_get(est, "cca_R2adj_env_total"); db <- com_get(est, "dbrda_R2adj_env_total")
  if (!is.na(cc) || !is.na(db)) out <- c(out, sprintf(
    "Other constrained ordinations of the same variables: CCA (chi-square distance) adjusted R2 = %s; db-RDA on %s dissimilarity adjusted R2 = %s. Different distances weight rare and common species differently, so these numbers need not match.",
    if (is.na(cc)) "not run" else pct(cc), com_dist_name(o$settings$distance), if (is.na(db)) "not run" else pct(db)))
  g <- data$group
  if (!is.null(g)) {
    pn <- paste0("R2adj_", g, "_permanova"); pp <- com_get(est, pn, "p_value")
    gr <- factor(data$env[[g]]); tb <- table(gr); bal <- length(unique(tb)) == 1
    s <- sprintf("%s groups (%s) %s in composition: group membership accounted for %s of the variation in %s dissimilarities (PERMANOVA %s, %s).",
                 lv(g), paste(sprintf("%s: %d sites", gsub("_", " ", names(tb)), as.integer(tb)), collapse = ", "),
                 if (!is.na(pp) && pp < 0.05) "differed" else "did not differ detectably",
                 ci(pn), com_dist_name(o$settings$distance), adjtxt, fmt_p(pp))
    dp <- com_get(est, paste0("dispersion_ratio_", g), "p_value")
    dr <- com_get(est, paste0("dispersion_ratio_", g))
    if (!is.na(dp)) {
      md <- if (!is.null(o$betadisper)) tapply(o$betadisper$distances, gr, mean) else NULL
      s <- paste(s, sprintf("Multivariate dispersion (mean distance to group centroid: %s) %s (betadisper permutation test %s; largest/smallest = %.2f).",
                            paste(sprintf("%s %.3f", gsub("_", " ", names(md)), md), collapse = ", "),
                            if (dp < 0.05) "differed among groups" else "did not differ detectably among groups", fmt_p(dp), dr))
      if (!is.na(pp) && pp < 0.05 && dp < 0.05)
        s <- paste(s, sprintf("Because dispersion differs, the PERMANOVA result cannot be read as a shift in average composition alone: PERMANOVA tests location and dispersion together, and %s. Inspect the ordination: if group centroids overlap but one group is more spread out, the groups differ in variability (heterogeneity among sites), not in typical species composition.",
                              if (bal) "with equal group sizes it is fairly robust, but part of the signal may be dispersion"
                              else "with unequal group sizes it can give false positives when the smaller group is the more variable one (Anderson & Walsh 2013)"))
      else if (!is.na(pp) && pp < 0.05) s <- paste(s, "Because PERMDISP found no dispersion difference, the PERMANOVA result is most consistent with a shift in average species composition (location). PERMDISP has limited power for modest differences in spread, though, so a smaller contribution from dispersion cannot be ruled out.")
      else if (!is.na(dp) && dp < 0.05) s <- paste(s, "The groups differ in how variable their communities are rather than in average composition.")
    }
    out <- c(out, s)
    pw <- o$permanova_pairwise
    if (!is.null(pw)) out <- c(out, paste0("Pairwise PERMANOVA (Holm-adjusted): ",
      paste(sprintf("%s R2 = %.3f, %s", gsub("_", " ", pw$pair), pw$R2, fmt_p(pw$p_holm)), collapse = "; "), "."))
    pe <- o$permanova_env
    if (!is.null(pe) && "g" %in% rownames(pe)) out <- c(out, sprintf(
      "Accounting for the environmental variables, %s still explained %s of the dissimilarity (marginal PERMANOVA, unadjusted R2, %s; for comparison, its unadjusted R2 without the environmental variables is %s).",
      tolower(lv(g)), pct(pe["g", "R2"]), fmt_p(pe["g", "Pr(>F)"]),
      if (!is.null(o$permanova)) pct(o$permanova[1, "R2"]) else "not available"))
    ind <- o$indicators
    if (!is.null(ind)) {
      sig <- ind[ind$p_fdr < 0.05, ]
      out <- c(out, if (nrow(sig)) paste0(sprintf("%d species %s of one %s after false-discovery-rate correction (IndVal.g; A = specificity, the share of the species' group-size-equalised mean abundance found in that group; B = fidelity, the share of that group's sites where it occurs): ",
                                                   nrow(sig), if (nrow(sig) == 1) "was an indicator" else "were indicators", tolower(lv(g))),
                                           paste(sprintf("%s for %s (IndVal %.2f, A %.2f, B %.2f, %s)", vapply(sig$species, lv, ""), gsub("_", " ", sig$group), sig$indval, sig$A, sig$B, fmt_p(sig$p_fdr)), collapse = "; "), ".")
             else sprintf("No species was a significant indicator of any %s after false-discovery-rate correction (%d had uncorrected p < 0.05, about %.1f expected by chance alone).",
                          tolower(lv(g)), sum(ind$p < 0.05), 0.05 * nrow(ind)))
    }
  }
  out <- c(out, sprintf(
    "Beta diversity (presence/absence, Baselga 2010): two sites differed on average by Sorensen dissimilarity %.2f (95%% CI %.2f to %.2f). Of this, %s (95%% CI %s to %s) was turnover (species replaced by other species) and the rest nestedness (species-poor sites holding subsets of richer sites' species). %s",
    com_get(est, "beta_sor_pair_mean"), com_get(est, "beta_sor_pair_mean", "lower"), com_get(est, "beta_sor_pair_mean", "upper"),
    pct(com_get(est, "turnover_share"), 0), pct(com_get(est, "turnover_share", "lower"), 0), pct(com_get(est, "turnover_share", "upper"), 0),
    if (com_get(est, "turnover_share") > 0.5) "Composition changes mainly by replacement: different sites hold different species, so no single species-rich site contains most of the others' species."
    else "Composition changes mainly by loss of species, so the poorer sites are largely subsets of the richer ones; the richest sites hold most species."))
  bm <- com_get(est, "beta_SOR_multisite")
  if (!is.na(bm)) out <- c(out, sprintf("Multiple-site Sorensen dissimilarity across all %d sites was %.2f (turnover %.2f, nestedness %.2f). Multiple-site values grow with the number of sites, so compare them only between data sets with the same number of sites (betapart::beta.sample).",
                                       nrow(data$comm), bm, com_get(est, "beta_SIM_multisite"), com_get(est, "beta_SNE_multisite")))
  stv <- com_get(est, "nmds_stress")
  if (!is.na(stv)) out <- c(out, sprintf("The 2-D NMDS has stress %.3f: %s (Clarke 1993: < 0.05 excellent, < 0.1 good, < 0.2 usable, > 0.2 close to arbitrary).", stv,
    if (stv < 0.05) "an excellent representation" else if (stv < 0.1) "a good representation" else if (stv < 0.2) "usable, but do not over-read small distances between points" else "a poor representation; use the constrained ordination and the tests, not the NMDS picture"))
  out <- c(out, "Dissimilarities and beta diversity describe the communities as sampled. Rare or hard-to-detect species are missed at some sites where they occur, which inflates apparent turnover; with unequal sampling effort among sites these numbers mix ecology and effort.")
  out
}

com_methods_text <- function(fit, data, settings) {
  st <- fit$object$settings
  sprintf(paste(
    "We analysed a site-by-species matrix of %d sites and %d species (%s).",
    "Species data were %s-transformed (Legendre & Gallagher 2001) and related to %s by redundancy analysis (RDA); variation was partitioned into unique and shared fractions using adjusted R2 (Borcard et al. 1992; Peres-Neto et al. 2006), with 95%% confidence intervals from a delete-one-site jackknife and significance of each testable fraction from %d permutations (vegan; Oksanen et al.).",
    "%s",
    "Beta diversity was partitioned into turnover (Simpson) and nestedness-resultant components of Sorensen dissimilarity (Baselga 2010) for all site pairs and across all sites (betapart; Baselga & Orme 2012); intervals for pairwise means are delete-one-site jackknife intervals.",
    "We report non-metric multidimensional scaling (NMDS; %s dissimilarity, %d dimensions, up to %d random starts) with environmental vectors fitted by envfit, detrended correspondence analysis (DCA) axis length to judge linear vs unimodal responses, canonical correspondence analysis (CCA; ter Braak 1986) and distance-based RDA (Legendre & Anderson 1999).",
    "Species response curves are quasi-Poisson (or quasi-binomial) GLMs quadratic in each gradient (ter Braak & Looman 1986). Analyses were run with ecotools module 'community' in R."),
    nrow(data$comm), ncol(data$comm), if (isTRUE(fit$object$presence)) "presence/absence" else "abundance",
    st$transform, paste(vapply(data$env_vars, function(v) tolower(com_lab(v, NULL, data)), ""), collapse = ", "), st$n_perm,
    if (is.null(data$group)) "" else sprintf(
      "Differences in composition among %s groups were tested by PERMANOVA (Anderson 2001; adonis2, %s dissimilarity, %d permutations), reported with adjusted R2 and a jackknife interval, followed by Holm-adjusted pairwise tests; homogeneity of multivariate dispersion was tested with PERMDISP (Anderson 2006; betadisper, distances to group centroids, %d permutations) because PERMANOVA is sensitive to dispersion differences in unbalanced designs (Anderson & Walsh 2013). Indicator species were identified with the group-equalised indicator value IndVal.g (Dufrene & Legendre 1997; De Caceres & Legendre 2009; indicspecies) with Benjamini-Hochberg false-discovery-rate correction.",
      tolower(com_lab(data$group, NULL, data)), com_dist_name(st$distance), st$n_perm, st$n_perm),
    com_dist_name(st$distance), st$nmds_k, st$nmds_trymax)
}

# ---- figures ----------------------------------------------------------------

com_arrow_mul <- function(arrows, pts, fill = 0.75) {
  r <- apply(abs(pts), 2, max); a <- apply(abs(arrows), 2, max)
  fill * min(r / pmax(a, 1e-12))
}

com_site_df <- function(sc, data, opts) {
  d <- data.frame(x = sc[, 1], y = sc[, 2])
  if (!is.null(data$group)) {
    g <- factor(data$env[[data$group]])
    levels(g) <- vapply(levels(g), function(l) fig_label(l, l, opts), "")
    d$group <- g
  }
  d
}

com_ord_plot <- function(sites, arrows, species, xlab, ylab, data, opts, ellipse = TRUE) {
  has_g <- "group" %in% names(sites)
  p <- ggplot2::ggplot(sites, ggplot2::aes(x, y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3)
  if (has_g && ellipse && all(table(sites$group) >= 4))
    p <- p + ggplot2::stat_ellipse(ggplot2::aes(colour = group), type = "t", level = 0.68, linewidth = 0.4, linetype = 2, show.legend = FALSE)
  p <- p + if (has_g) ggplot2::geom_point(ggplot2::aes(colour = group, shape = group), size = 1.6, alpha = 0.85)
           else ggplot2::geom_point(colour = eco_pal$accent, size = 1.6, alpha = 0.85)
  has_sp <- !is.null(species) && nrow(species); has_ar <- !is.null(arrows) && nrow(arrows)
  if (has_sp)
    p <- p + ggplot2::geom_point(data = species, ggplot2::aes(x, y), inherit.aes = FALSE, shape = 3, size = 0.8, colour = eco_pal$muted)
  # Species labels are repelled; arrow labels sit just outside each arrow tip
  # (justified outward). Empty-label obstacles along each arrow and at its tip keep the
  # species labels off the arrows and their labels.
  if (has_sp) {
    obst <- if (has_ar) do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1, 1.1), function(f)
      data.frame(x = f * arrows$x, y = f * arrows$y, label = ""))) else NULL
    lab_df <- rbind(data.frame(x = species$x, y = species$y, label = species$label), obst)
    p <- p + ggrepel::geom_text_repel(data = lab_df, ggplot2::aes(x, y, label = label), inherit.aes = FALSE,
                                      size = 2.0, colour = eco_pal$muted, fontface = "italic", max.overlaps = 30,
                                      segment.size = 0.2, min.segment.length = 0.2, seed = 1, box.padding = 0.3)
  }
  if (has_ar) {
    # approximate label extent in data units, so outward labels are not clipped
    xr <- diff(range(c(sites$x, arrows$x, if (has_sp) species$x))); yr <- diff(range(c(sites$y, arrows$y, if (has_sp) species$y)))
    w <- nchar(arrows$label) * 0.02 * xr
    room <- data.frame(x = arrows$x * 1.04 + ifelse(arrows$x >= 0, 1, -1) * w,
                       y = arrows$y * 1.04 + ifelse(arrows$y >= 0, 1, -1) * 0.06 * yr)
    p <- p + ggplot2::geom_segment(data = arrows, ggplot2::aes(x = 0, y = 0, xend = x, yend = y), inherit.aes = FALSE,
                                   arrow = ggplot2::arrow(length = ggplot2::unit(1.6, "mm")), colour = eco_pal$ink, linewidth = 0.5) +
      ggplot2::geom_text(data = arrows, ggplot2::aes(x * 1.04, y * 1.04, label = label), inherit.aes = FALSE,
                         size = 2.4, colour = eco_pal$ink, fontface = "bold",
                         hjust = ifelse(abs(arrows$x) < 0.3 * abs(arrows$y), 0.5, ifelse(arrows$x >= 0, 0, 1)),
                         vjust = ifelse(abs(arrows$y) < 0.3 * abs(arrows$x), 0.5, ifelse(arrows$y >= 0, 0, 1))) +
      ggplot2::geom_blank(data = room, ggplot2::aes(x, y), inherit.aes = FALSE)  # room for the text
  }
  if (has_g) p <- p + scale_colour_eco(opts$palette %||% "okabe_ito")
  p <- p + ggplot2::coord_equal() + ggplot2::labs(x = xlab, y = ylab, colour = com_lab(data$group %||% "group", opts, data),
                                             shape = com_lab(data$group %||% "group", opts, data))
  attr(p, "fig_size") <- list(width = "double", aspect = 0.62)
  p
}

com_wrap <- function(w) function(x) vapply(strwrap(x, w, simplify = FALSE), paste, "", collapse = "\n")

com_top_species <- function(spsc, n) {
  if (is.null(spsc) || !nrow(spsc) || n <= 0) return(NULL)
  r <- sqrt(rowSums(spsc[, 1:2, drop = FALSE]^2))
  spsc[order(r, decreasing = TRUE)[seq_len(min(n, nrow(spsc)))], , drop = FALSE]
}

com_figures <- function(fit, data, opts = fig_opts()) {
  eco_require(c("ggplot2", "ggrepel"))
  o <- fit$object; est <- fit$estimates; figs <- list()
  ci <- opts$ci %||% 0.95
  n_sp <- if (is.null(opts$n_species_labels)) 8 else opts$n_species_labels
  lab <- function(k) com_lab(k, opts, data)
  # 1. NMDS ordination with envfit vectors
  if (!is.null(o$nmds)) {
    sc <- vegan::scores(o$nmds, display = "sites")
    sites <- com_site_df(sc, data, opts)
    arrows <- NULL
    if (!is.null(o$envfit) && !is.null(o$envfit$vectors)) {
      v <- o$envfit$vectors
      a <- sweep(v$arrows, 1, sqrt(v$r), "*")
      a <- a * com_arrow_mul(a, sc)
      arrows <- data.frame(x = a[, 1], y = a[, 2], label = vapply(rownames(a), lab, ""))
      arr_txt <- paste(sprintf("%s (r2 = %.2f, %s)", arrows$label, v$r, vapply(v$pvals, fmt_p, "")), collapse = "; ")
    } else arr_txt <- "none"
    spsc <- tryCatch(vegan::scores(o$nmds, display = "species"), error = function(e) NULL)
    spd <- com_top_species(spsc, n_sp)
    if (!is.null(spd)) spd <- data.frame(x = spd[, 1], y = spd[, 2], label = vapply(rownames(spd), lab, ""))
    figs$ordination <- with_caption(
      com_ord_plot(sites, arrows, spd, "NMDS 1", "NMDS 2", data, opts),
      sprintf("Non-metric multidimensional scaling of %s dissimilarities (stress %.3f). Points are sites%s; dashed ellipses enclose about 68%% of each group's sites. Arrows: environmental variables fitted to the ordination (envfit), length proportional to correlation: %s. Grey crosses and italic labels: weighted-average positions of the %d species furthest from the centre.",
              com_dist_name(o$settings$distance), o$nmds$stress, if (!is.null(data$group)) paste0(" coloured by ", tolower(lab(data$group))) else "",
              arr_txt, if (is.null(spd)) 0 else nrow(spd)))
    # Shepard diagram
    sh <- data.frame(diss = o$nmds$diss, dist = o$nmds$dist, dhat = o$nmds$dhat)
    sh <- sh[order(sh$diss), ]
    figs$nmds_shepard <- with_caption(
      ggplot2::ggplot(sh, ggplot2::aes(diss, dist)) +
        ggplot2::geom_point(size = 0.4, alpha = 0.25, colour = eco_pal$muted) +
        ggplot2::geom_step(ggplot2::aes(y = dhat), colour = eco_pal$truth, linewidth = 0.6) +
        ggplot2::labs(x = sprintf("Observed %s dissimilarity", com_dist_name(o$settings$distance)), y = "Distance in the NMDS ordination"),
      sprintf("Shepard diagram: each point is a pair of sites. The red step line is the monotone fit; scatter around it is what stress measures (stress %.3f).", o$nmds$stress))
  }
  # 2. RDA triplot
  if (!is.null(o$rda)) {
    r <- o$rda
    ax <- if (r$CCA$rank >= 2) c(1, 2) else c(1, 1)
    sc_all <- vegan::scores(r, choices = if (r$CCA$rank >= 2) 1:2 else 1:2, display = c("sites", "species", "bp"), scaling = 2)
    sc <- sc_all$sites
    sites <- com_site_df(sc, data, opts)
    bp <- sc_all$biplot
    bp <- bp * com_arrow_mul(bp, sc)
    arrows <- data.frame(x = bp[, 1], y = bp[, 2], label = vapply(rownames(bp), lab, ""))
    spd <- com_top_species(sc_all$species, n_sp)
    if (!is.null(spd)) spd <- data.frame(x = spd[, 1], y = spd[, 2], label = vapply(rownames(spd), lab, ""))
    eig <- r$CCA$eig; tot <- r$tot.chi
    axlab <- function(k) {
      nm <- colnames(sc)[k]
      ev <- if (nm %in% names(eig)) eig[[nm]] else r$CA$eig[[nm]]
      sprintf("%s (%.1f%% of total variation)", nm, 100 * ev / tot)
    }
    figs$rda_triplot <- with_caption(
      com_ord_plot(sites, arrows, spd, axlab(1), axlab(2), data, opts),
      sprintf("Redundancy analysis of %s-transformed species data constrained by %s (scaling 2: angles between arrows reflect correlations). Points are sites; arrows are the environmental variables; italic labels are the %d species best represented in this plane. Constrained axes together explain adjusted R2 = %s of the variation.",
              o$settings$transform, paste(vapply(data$env_vars, lab, ""), collapse = ", "), if (is.null(spd)) 0 else nrow(spd),
              pct(com_get(est, "R2adj_env_total"))))
  }
  # 3. Variance partition
  part <- names(o$partition)
  if (length(part)) {
    rows <- if (length(part) > 1) c(paste0("R2adj_unique_", part), "R2adj_shared", "R2adj_env_total") else "R2adj_env_total"
    vp <- est[match(rows, est$param), c("param", "estimate", "lower", "upper")]
    vp$label <- c(if (length(part) > 1) vapply(part, function(k) paste(lab(k), "(unique)"), ""),
                  if (length(part) > 1) "Shared", "All variables")
    vp <- rbind(vp, data.frame(param = "residual", estimate = 1 - com_get(est, "R2adj_env_total"), lower = NA, upper = NA,
                               label = "Unexplained"))
    vp$label <- factor(vp$label, levels = rev(vp$label))
    vp$kind <- ifelse(vp$param == "R2adj_env_total", "total", ifelse(vp$param == "residual", "residual", "fraction"))
    pvp <- ggplot2::ggplot(vp, ggplot2::aes(estimate, label)) +
        ggplot2::geom_col(ggplot2::aes(fill = kind), width = 0.65, show.legend = FALSE) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = 0.25, linewidth = 0.4, colour = eco_pal$ink, na.rm = TRUE) +
        ggplot2::geom_vline(xintercept = 0, colour = eco_pal$ink, linewidth = 0.3) +
        ggplot2::scale_fill_manual(values = c(fraction = "#56B4E9", total = "#0072B2", residual = "grey80")) +
        ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
        ggplot2::labs(x = "Share of variation (adjusted R\u00b2)", y = NULL) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1))
    attr(pvp, "fig_size") <- list(width = "onehalf", aspect = 0.55)
    figs$variance_partition <- with_caption(pvp,
      sprintf("Variation partitioning (RDA, %s transformation). Bars: adjusted R2 of each fraction; error bars: %d%% jackknife confidence intervals. 'Unique' = explained by that set after removing the others; 'Shared' = jointly explained and not attributable to either set (can be negative). Unexplained = 1 - adjusted R2 of all variables.",
              o$settings$transform, round(100 * ci)))
  }
  # 4. Effect sizes (all adjusted R2 with CIs)
  es <- est[grepl("^R2adj_", est$param), c("param", "estimate", "lower", "upper")]
  if (nrow(es)) {
    nice <- function(p) {
      k <- sub("^R2adj_", "", p)
      if (k == "env_total") return("All environmental variables")
      if (k == "shared") return("Shared (environment)")
      if (grepl("^unique_", k)) return(paste(lab(sub("^unique_", "", k)), "(unique)"))
      if (grepl("_permanova$", k)) return(paste(lab(sub("_permanova$", "", k)), "(PERMANOVA)"))
      lab(k)
    }
    es$param <- vapply(es$param, function(p) fig_label(p, nice(p), opts), "")
    pe <- fig_coef(es, null_value = 0, xlab = sprintf("Adjusted R\u00b2 (%d%% CI)", round(100 * ci)), opts = NULL) +
      ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"))
    attr(pe, "fig_size") <- list(width = "onehalf", aspect = 0.6)
    figs$effect_sizes <- with_caption(pe,
      "Share of variation in community composition associated with each environmental set (RDA) and with the grouping factor (PERMANOVA). Dashed line: no association. Points: adjusted R2; lines: jackknife confidence intervals.")
  }
  # 5. Beta components
  bp <- o$beta_pairs
  if (!is.null(bp)) {
    ut <- upper.tri(bp$sor)
    comp_lab <- c(sor = fig_label("beta_sor", "Total (Sorensen)", opts), sim = fig_label("beta_sim", "Turnover (Simpson)", opts),
                  sne = fig_label("beta_sne", "Nestedness", opts))
    pd <- do.call(rbind, lapply(names(comp_lab), function(k) data.frame(component = comp_lab[[k]], value = bp[[k]][ut])))
    pd$component <- factor(pd$component, levels = comp_lab)
    mn <- data.frame(component = factor(comp_lab, levels = comp_lab),
                     estimate = vapply(names(comp_lab), function(k) com_get(est, paste0("beta_", k, "_pair_mean")), 0),
                     lower = vapply(names(comp_lab), function(k) com_get(est, paste0("beta_", k, "_pair_mean"), "lower"), 0),
                     upper = vapply(names(comp_lab), function(k) com_get(est, paste0("beta_", k, "_pair_mean"), "upper"), 0))
    figs$beta_components <- with_caption(
      ggplot2::ggplot(pd, ggplot2::aes(component, value)) +
        ggplot2::geom_violin(fill = "grey90", colour = "grey70", linewidth = 0.3, scale = "width", bw = 0.05) +
        ggplot2::geom_pointrange(data = mn, ggplot2::aes(y = estimate, ymin = lower, ymax = upper), colour = eco_pal$ink, size = 0.3) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::scale_x_discrete(labels = com_wrap(12)) +
        ggplot2::labs(x = NULL, y = "Pairwise dissimilarity"),
      sprintf("Baselga (2010) partition of pairwise Sorensen dissimilarity into turnover (Simpson) and nestedness-resultant components. Grey: distribution over all %d site pairs. Points and bars: mean over pairs with %d%% jackknife CI. Turnover share = %s.",
              sum(ut), round(100 * ci), pct(com_get(est, "turnover_share"), 0)))
    # distance decay along the first numeric gradient
    nv <- data$env_vars[vapply(data$env_vars, function(v) is.numeric(data$env[[v]]), TRUE)]
    if (length(nv)) {
      v <- nv[1]; x <- data$env[[v]]
      dx <- abs(outer(x, x, "-"))[ut]
      dd <- rbind(data.frame(dx = dx, value = bp$sim[ut], component = comp_lab[["sim"]]),
                  data.frame(dx = dx, value = bp$sne[ut], component = comp_lab[["sne"]]))
      pdd <- ggplot2::ggplot(dd, ggplot2::aes(dx, value, colour = component, linetype = component)) +
          ggplot2::geom_point(size = 0.3, alpha = 0.15) +
          ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8, ggplot2::aes(fill = component), alpha = 0.2) +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          scale_colour_eco(opts$palette %||% "okabe_ito") + scale_fill_eco(opts$palette %||% "okabe_ito") +
          ggplot2::labs(x = paste("Difference in", tolower(lab(v))), y = "Pairwise dissimilarity", colour = NULL, linetype = NULL, fill = NULL)
      attr(pdd, "fig_size") <- list(width = "onehalf", aspect = 0.6)
      figs$distance_decay <- with_caption(pdd,
        sprintf("Turnover and nestedness components of dissimilarity for every pair of sites against their difference in %s. Lines: loess smooths with 95%% confidence bands. Turnover rising with environmental distance indicates replacement of species along the gradient.",
                tolower(lab(v))))
    }
  }
  # 6. Dispersion
  if (!is.null(o$betadisper)) {
    bd <- o$betadisper
    g <- factor(data$env[[data$group]]); levels(g) <- vapply(levels(g), function(l) fig_label(l, l, opts), "")
    dd <- data.frame(group = g, dist = bd$distances)
    z <- stats::qt(1 - (1 - ci) / 2, pmax(table(g) - 1, 1))
    ms <- do.call(rbind, lapply(levels(g), function(l) { v <- dd$dist[dd$group == l]; se <- stats::sd(v) / sqrt(length(v))
      data.frame(group = l, estimate = mean(v), lower = mean(v) - z[[l]] * se, upper = mean(v) + z[[l]] * se) }))
    ms$group <- factor(ms$group, levels = levels(g))
    pp <- o$betadisper_test$tab[1, "Pr(>F)"]
    figs$dispersion <- with_caption(
      ggplot2::ggplot(dd, ggplot2::aes(group, dist, colour = group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey60", linewidth = 0.3) +
        ggplot2::geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.6) +
        ggplot2::geom_pointrange(data = ms, ggplot2::aes(y = estimate, ymin = lower, ymax = upper), colour = eco_pal$ink,
                                 size = 0.3, position = ggplot2::position_nudge(x = 0.32)) +
        scale_colour_eco(opts$palette %||% "okabe_ito") + ggplot2::guides(colour = "none") +
        ggplot2::scale_x_discrete(labels = com_wrap(12)) +
        ggplot2::labs(x = lab(data$group), y = "Distance to group centroid"),
      sprintf("Multivariate dispersion (PERMDISP, betadisper on %s dissimilarity): distance of each site (points) to its group's centroid. Boxes: quartiles; black points and bars: group means with %d%% t confidence intervals. Permutation test of equal dispersion: %s. PERMANOVA cannot distinguish location from dispersion differences, so read it together with this figure.",
              com_dist_name(o$settings$distance), round(100 * ci), fmt_p(pp)))
  }
  # 7. Indicator species
  if (!is.null(o$indicators) && nrow(o$indicators)) {
    ind <- o$indicators
    sel <- ind[ind$p < 0.05, ]
    if (!nrow(sel)) sel <- utils::head(ind, 10)
    sel <- utils::head(sel[order(sel$group, -sel$indval), ], 25)
    sel$label <- vapply(sel$species, lab, "")
    sel$label <- factor(sel$label, levels = rev(unique(sel$label)))
    sel$group_l <- factor(vapply(sel$group, function(l) fig_label(l, l, opts), ""),
                          levels = vapply(levels(factor(data$env[[data$group]])), function(l) fig_label(l, l, opts), ""))
    sel$sig <- ifelse(sel$p_fdr < 0.05, "FDR < 0.05", "p < 0.05 (uncorrected)")
    pin <- ggplot2::ggplot(sel, ggplot2::aes(indval, label, colour = group_l, shape = sig)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = indval, yend = label), linewidth = 0.4) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_shape_manual(values = c("FDR < 0.05" = 16, "p < 0.05 (uncorrected)" = 1), drop = FALSE) +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
        scale_colour_eco(opts$palette %||% "okabe_ito", drop = FALSE) +
        ggplot2::labs(x = "Indicator value (IndVal.g, square root)", y = NULL, colour = lab(data$group), shape = NULL) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
    attr(pin, "fig_size") <- list(width = "double", height = 30 + 6 * nrow(sel))
    figs$indicator_species <- with_caption(pin,
      sprintf("Indicator species analysis (IndVal.g, %d permutations). Each species is shown for the %s it is most associated with; IndVal combines specificity (share of the species' abundance in that group) and fidelity (share of that group's sites occupied). Filled: significant after Benjamini-Hochberg correction; open: uncorrected p < 0.05 only.%s",
              o$settings$n_perm, tolower(lab(data$group)), if (all(ind$p >= 0.05)) " No species had p < 0.05; the 10 highest values are shown." else ""))
  }
  # 8. Species response curves along each gradient (response scale)
  if (!is.null(o$curves)) for (v in names(o$curves)) {
    cv <- o$curves[[v]]$pred
    cv$species_l <- factor(vapply(cv$species, lab, ""), levels = unique(vapply(cv$species, lab, "")))
    obs <- data.frame(x = data$env[[v]])
    p <- ggplot2::ggplot(cv, ggplot2::aes(x, estimate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.2) +
      ggplot2::geom_line(colour = eco_pal$band, linewidth = 0.7) +
      ggplot2::geom_rug(data = obs, ggplot2::aes(x = x), inherit.aes = FALSE, alpha = 0.3, sides = "b", length = ggplot2::unit(0.03, "npc")) +
      ggplot2::facet_wrap(~species_l, scales = "free_y", ncol = 3) +
      ggplot2::labs(x = lab(v), y = if (isTRUE(o$presence)) "Probability of occurrence" else "Expected count per site") +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic"))
    if (isTRUE(o$presence)) p <- p + ggplot2::scale_y_continuous(limits = c(0, 1))
    if (isTRUE(opts$show_data)) {
      od <- do.call(rbind, lapply(unique(cv$species), function(sp) data.frame(x = data$env[[v]], estimate = if (isTRUE(o$presence)) (data$comm[, sp] > 0) * 1 else data$comm[, sp],
                                                                         species_l = factor(lab(sp), levels = levels(cv$species_l)))))
      p <- p + ggplot2::geom_point(data = od, size = 0.5, alpha = 0.3, colour = eco_pal$muted)
    }
    attr(p, "fig_size") <- list(width = "double", aspect = 0.55)
    figs[[paste0("species_response_", v)]] <- with_caption(p, sprintf(
      "Response curves of the %d species whose occurrence or abundance is most strongly associated with %s (quadratic %s GLM; other numeric gradients held at their means). Lines: fitted %s; bands: %d%% confidence intervals; points: observed sites; rug: sampled values of %s.",
      length(unique(cv$species)), tolower(lab(v)), if (isTRUE(o$presence)) "quasi-binomial" else "quasi-Poisson",
      if (isTRUE(o$presence)) "probability of occurrence" else "expected count", round(100 * ci), tolower(lab(v))))
  }
  figs
}
