# Interpretation, methods text and publication figures for popgen-structure.

ps_row <- function(est, p) est[est$param == p, , drop = FALSE]
ps_fmt <- function(r, d = 3) sprintf(paste0("%.", d, "f (95%% CI %.", d, "f to %.", d, "f)"), r$estimate, r$lower, r$upper)

ps_wright_label <- function(f) {
  if (f < 0.05) "little" else if (f < 0.15) "moderate" else if (f < 0.25) "great" else "very great"
}

ps_interpret <- function(fit, data) {
  e <- fit$estimates; o <- fit$object; pops <- o$comps$pops
  out <- character()
  fst <- ps_row(e, "Fst_global")
  excl0 <- fst$lower > 0
  out <- c(out, sprintf(
    "Overall differentiation among the %d sampled populations is FST = %s (Weir & Cockerham theta). %s",
    length(pops), ps_fmt(fst),
    if (excl0) sprintf("The interval excludes 0, so the samples do not come from one randomly mating population. On Wright's descriptive scale this is %s differentiation: about %.1f%% of the gene diversity lies between populations.",
                       ps_wright_label(fst$estimate), 100 * fst$estimate)
    else "The interval includes 0: these data cannot distinguish the sampled sites from a single randomly mating population (lack of evidence for structure is not proof of connectivity; see the power analysis for how many loci and individuals would detect weak structure)."))
  nm <- ps_row(e, "Nm_island_model")
  if (excl0 && is.finite(nm$estimate))
    out <- c(out, sprintf(
      "If the populations behaved like Wright's idealised island model at equilibrium, this FST would correspond to about %.1f effective migrants per generation (range %.1f to %.1f across the FST interval). Treat this only as an order of magnitude: real populations violate the island model (unequal sizes, non-equilibrium after recent fragmentation, isolation by distance), and FST is a non-linear function of Nm, so small FST errors become large Nm errors (Whitlock & McCauley 1999). Mills & Allendorf (1996) treat one migrant per generation as a minimum and suggest 1 to 10 migrants per generation as a general guideline for limiting loss of diversity by drift.",
      nm$estimate, nm$lower, nm$upper))
  pw <- e[grepl("^Fst_.*_vs_", e$param), ]
  if (nrow(pw) >= 2) {
    hi <- pw[which.max(pw$estimate), ]; lo <- pw[which.min(pw$estimate), ]
    lab <- function(p) gsub("_vs_", " and ", sub("^Fst_", "", p))
    nz <- sum(pw$lower > 0)
    out <- c(out, sprintf(
      "Pairwise FST ranges from %s between %s to %s between %s; %d of %d pairs have intervals excluding 0. Highly differentiated pairs exchange few migrants relative to drift; they are the places where connectivity is weakest (or where one population is small and drifting fast).",
      ps_fmt(lo), lab(lo$param), ps_fmt(hi), lab(hi$param), nz, nrow(pw)))
  }
  he <- e[grepl("^He_", e$param) & e$param != "He_mean", ]; ar <- e[grepl("^AR_", e$param), ]
  hm <- ps_row(e, "He_mean")
  low <- he[which.min(he$estimate), ]; high <- he[which.max(he$estimate), ]
  arl <- ar[which.min(ar$estimate), ]
  flag <- he$param[he$upper < hm$estimate]
  out <- c(out, sprintf(
    "Genetic diversity: mean expected heterozygosity is %s. It is lowest in %s (%s) and highest in %s (%.3f). Allelic richness (rarefied to %d gene copies so populations with different sample sizes are comparable) is lowest in %s (%.2f alleles per locus). %s",
    ps_fmt(hm), sub("^He_", "", low$param), ps_fmt(low), sub("^He_", "", high$param), high$estimate,
    o$comps$rarefaction_genes, sub("^AR_", "", arl$param), arl$estimate,
    if (length(flag)) sprintf("Population(s) %s have He significantly below the overall mean (their interval lies below it); small, isolated populations lose diversity by drift and are the first candidates to monitor for inbreeding depression and to consider for genetic rescue (augmented gene flow). Frankham et al. (2011) found a low risk of outbreeding depression when source and recipient are the same species, have the same karyotype, have been isolated for fewer than about 500 years and live in similar environments; a large FST alone is not evidence of local adaptation.",
                              paste(sub("^He_", "", flag), collapse = ", "))
    else "No population has He clearly below the overall mean."))
  fis <- ps_row(e, "Fis_global")
  fp <- e[grepl("^Fis_", e$param) & e$param != "Fis_global", ]
  pos <- sub("^Fis_", "", fp$param[fp$lower > 0]); neg <- sub("^Fis_", "", fp$param[fp$upper < 0])
  out <- c(out, sprintf(
    "Within populations, FIS = %s. %s%s%s%s",
    ps_fmt(fis),
    if (fis$lower > 0) "Individuals are more homozygous than expected under random mating. Possible causes are inbreeding (mating among relatives or selfing), null alleles (microsatellites), or sampling more than one breeding group as one population (Wahlund effect); check whether the excess is spread over all loci (inbreeding) or concentrated in a few (null alleles, scoring). "
    else if (fis$upper < 0) "Individuals are slightly more heterozygous than expected, as happens in very small populations or with avoidance of mating between relatives. "
    else "Genotype proportions are consistent with random mating within populations. ",
    if (length(pos)) sprintf("Heterozygote deficit (FIS interval above 0) in: %s. ", paste(pos, collapse = ", ")) else "",
    if (length(neg)) sprintf("Heterozygote excess in: %s. ", paste(neg, collapse = ", ")) else "",
    if (length(c(pos, neg))) sprintf("With %d populations each tested separately, about %.1f interval(s) would exclude 0 by chance alone even under random mating, so treat single-population FIS signals as tentative unless they are large or supported by the Hardy-Weinberg tests.", nrow(fp), 0.05 * nrow(fp)) else ""))
  if (!is.null(o$hwe)) {
    h <- o$hwe
    out <- c(out, sprintf("Hardy-Weinberg tests: %d of %d population x locus tests deviate after false-discovery-rate correction (q < 0.05).",
                          sum(h$q < 0.05, na.rm = TRUE), sum(!is.na(h$q))))
  }
  if (!is.null(o$amova)) {
    am <- o$amova
    out <- c(out, sprintf(
      "AMOVA (allele-distance based, poppr): %.1f%% of molecular variance lies between populations (PhiST = %.3f, permutation P = %.3f). PhiST uses allele-identity distances among individuals and is not identical to the Weir & Cockerham FST above; both measure the same kind of structure.",
      am$table[1, 2], am$phi[1, 1], am$p[1]))
  }
  sl <- ps_row(e, "ibd_slope")
  if (nrow(sl)) {
    lg <- identical(o$ibd_distance, "log")
    mt <- o$mantel
    out <- c(out, sprintf(
      "Isolation by distance: FST/(1-FST) increases by %s per %s (regression over population pairs)%s. %s",
      ps_fmt(sl, 4), if (lg) "unit of log distance (i.e. per 2.7-fold increase in distance)" else "km",
      if (!is.null(mt)) sprintf("; Mantel r = %.2f, permutation P = %.3f", mt$r, mt$p) else "",
      if (sl$lower > 0 && !is.null(mt) && mt$p >= 0.05) ps_ibd_mixed(e, o, mt)
      else if (sl$lower > 0) paste0("Differentiation grows with distance, consistent with gene flow that is mostly among nearby populations rather than across the whole study area; nearby populations are more connected than distant ones, so maintaining stepping-stone habitat between sites matters for connectivity. ",
                               if (lg) sprintf("Under Rousset's (1997) two-dimensional model the slope estimates 1/(4 pi D sigma^2), giving a neighbourhood size of roughly %.0f (%.0f to %.0f) breeding individuals, where D is effective density and sigma the parent-offspring dispersal distance.",
                                               1 / sl$estimate, 1 / sl$upper, 1 / sl$lower)
                               else sprintf("Under Rousset's (1997) one-dimensional model the slope estimates 1/(4 D sigma^2) per km (D = effective linear density, sigma = dispersal distance), i.e. 4 D sigma^2 ~ %.0f km.", 1 / sl$estimate))
      else "The interval includes 0: no clear increase of differentiation with distance in these data (with few populations the test has little power)."))
  }
  if (!is.null(o$dapc)) {
    dp <- o$dapc
    tabk <- table(data$pop, dp$grp)
    out <- c(out, sprintf(
      "Clustering without using population labels (k-means on principal components, number of clusters K chosen by BIC) found K = %d genetic cluster%s for %d sampled populations. %s Clusters are a summary of allele-frequency differences, not proof of distinct management units; with weak FST or isolation by distance, clustering methods can merge or split populations arbitrarily.",
      dp$K, if (dp$K > 1) "s" else "", length(pops),
      if (dp$K == 1) "The data are best described as one genetic group." else if (dp$K < length(pops)) "Some sampled populations share a cluster: they are genetically similar enough to be grouped (see the assignment plot)." else "Clusters correspond roughly to sampled populations (see the assignment plot)."))
  }
  pa <- e[grepl("^PrivateAlleles_", e$param), ]
  if (sum(pa$estimate) > 0) {
    top <- pa[which.max(pa$estimate), ]
    out <- c(out, sprintf("Private alleles (found in only one sampled population): %d in total, most in %s (%d). Private alleles mark unique variation worth conserving, but their count rises with sample size and they are often rare.",
                          as.integer(sum(pa$estimate)), sub("^PrivateAlleles_", "", top$param), as.integer(top$estimate)))
  }
  out <- c(out, "All statements describe associations in the genetic data. Gene flow, drift and population history all shape FST; these numbers do not by themselves identify what caused the structure.")
  out
}

ps_methods_text <- function(fit, data, settings) {
  o <- fit$object
  sprintf(paste(
    "Genotypes of %d individuals from %d populations at %d %s loci were analysed with the ecotools popgen-structure module (R %s; hierfstat-equivalent estimators verified against hierfstat %s).",
    "Observed heterozygosity (Ho) and Nei's (1987) unbiased gene diversity (He) were computed per population and locus and averaged over loci; allelic richness was rarefied to %d gene copies (El Mousadik & Petit 1996).",
    "FST and FIS were estimated with Weir & Cockerham's (1984) theta and f as ratios of variance components summed over alleles and loci, globally and for every pair of populations.",
    "95%% confidence intervals were estimate +/- t(L-1) x SE (L = %d loci), with SE^2 the delete-one-locus jackknife variance (Weir 1996) plus the covariance among loci caused by sharing individuals, estimated by a delete-one-individual jackknife stratified by population (joint minus locus-by-locus deletion variance).",
    "%s%s%s%s"),
    length(data$ind), nlevels(data$pop), length(data$loci), if (identical(data$marker, "snp")) "SNP" else "microsatellite",
    paste(R.version$major, R.version$minor, sep = "."), as.character(utils::packageVersion("hierfstat")),
    o$comps$rarefaction_genes, length(data$loci),
    if (!is.null(o$hwe)) sprintf("Departures from Hardy-Weinberg proportions were tested per population and locus with exact tests (biallelic; Wigginton et al. 2005) or Monte Carlo conditional exact tests (%d permutations; Guo & Thompson 1992), with Benjamini-Hochberg correction within populations. ", settings$hwe_B %||% 1000) else "",
    if (!is.null(o$dist)) sprintf("Isolation by distance was assessed by regressing pairwise FST/(1-FST) on %s geographic distance (Rousset 1997), with a confidence interval for the slope computed as above%s. ",
                                  if (identical(o$ibd_distance, "log")) "the natural logarithm of" else "",
                                  if (!is.null(o$mantel)) sprintf(" and a Mantel test (%d permutations; vegan)", o$mantel$perm) else "") else "",
    if (!is.null(o$pca) || !is.null(o$dapc)) "Population structure was summarised by principal component analysis of centred allele frequencies (adegenet/ade4; Jombart 2008); the number of genetic clusters was chosen by BIC from k-means clustering on principal components, and memberships were estimated by discriminant analysis of principal components (DAPC; Jombart et al. 2010) retaining K-1 principal components (Thia 2023). " else "",
    if (!is.null(o$amova)) "Hierarchical partitioning of molecular variance used AMOVA (Excoffier et al. 1992; poppr). " else "")
}

# ---- figures ------------------------------------------------------------------
#
# Figure ids (docs/FIGURES.md): diversity, fst_matrix, ibd, pca, bic_k,
# dapc_assignment, sample_map. Every label goes through fig_label(); population
# names can be renamed with `labels: {<POP>: "Name"}`. Intervals are redrawn at
# opts$ci from the estimates' SE (estimate +/- t(L-1) x SE, the engine's own
# interval formula), so `figures: {ci: 0.9}` changes the figures without refitting.
# Colour follows the population: the same population has the same colour in
# pca and dapc_assignment (a cluster made up mainly of one population takes
# its colour; shared clusters are grey); shape is the second encoding for
# greyscale. sample_map uses viridis fill for He (a magnitude), not population colour.

# Interval at level `ci` from the SE column (rows without SE keep NA).
ps_ci_rows <- function(e, ci, df) {
  tq <- stats::qt(1 - (1 - ci) / 2, max(df, 1))
  e$lower <- e$estimate - tq * e$sd; e$upper <- e$estimate + tq * e$sd
  e
}

# Fixed colour per population (entity), in the order the populations are listed.
ps_pop_cols <- function(pops, opts) {
  n <- length(pops)
  pal <- opts$palette %||% "okabe_ito"
  if (isTRUE(opts$greyscale) || identical(pal, "grey"))
    return(stats::setNames(grDevices::grey.colors(n, start = 0.05, end = 0.8), pops))
  p <- if (length(pal) > 1) pal else eco_pal[[pal]]
  if (n <= length(p)) stats::setNames(p[seq_len(n)], pops)
  else stats::setNames(grDevices::hcl.colors(n, "viridis"), pops)
}
ps_pop_shapes <- function(pops) {
  sh <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6)
  stats::setNames(if (length(pops) <= length(sh)) sh[seq_along(pops)] else rep(16, length(pops)), pops)
}

# Kept for backward compatibility (scale by name for any vector of groups).
ps_pop_scale <- function(pops, opts, aes = "colour") {
  vals <- ps_pop_cols(pops, opts)
  if (aes == "colour") ggplot2::scale_colour_manual(values = vals, breaks = pops)
  else ggplot2::scale_fill_manual(values = vals, breaks = pops)
}

# Cluster colours follow the population colours used in the PCA: a cluster in
# which exactly one population has mean membership > 0.5 (and is that
# population's best cluster) takes that population's colour; a cluster shared
# by several populations is drawn in neutral grey (darker greys for further
# shared clusters), so it is not mistaken for any one population.
ps_cluster_cols <- function(post, pop_lab, popcols) {
  cl <- colnames(post)
  mm <- apply(post, 2, function(p) tapply(p, pop_lab, mean))
  mm <- matrix(mm, nrow = nlevels(pop_lab), dimnames = list(levels(pop_lab), cl))
  out <- stats::setNames(rep(NA_character_, length(cl)), cl)
  for (j in seq_along(cl)) {
    maj <- rownames(mm)[mm[, j] > 0.5]
    if (length(maj) == 1 && which.max(mm[maj, ]) == j) out[j] <- popcols[[maj]]
  }
  greys <- c("#BDBDBD", "#8C8C8C", "#5E5E5E", "#D9D9D9", "#3B3B3B")
  spare <- c(greys, setdiff(eco_pal$okabe_ito, out), grDevices::hcl.colors(length(cl), "viridis"))
  out[is.na(out)] <- spare[seq_len(sum(is.na(out)))]
  out
}

ps_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  o <- fit$object; pops <- o$comps$pops
  L <- length(data$loci)
  ci_lvl <- opts$ci %||% o$settings$ci %||% 0.95
  ci <- round(100 * ci_lvl)
  e <- ps_ci_rows(fit$estimates, ci_lvl, L - 1)
  plab <- vapply(pops, function(p) fig_label(p, p, opts), "")
  names(plab) <- pops
  lab_pop <- fig_label("population", "Population", opts)
  popcols <- ps_pop_cols(unname(plab), opts); popshape <- ps_pop_shapes(unname(plab))
  figs <- list()

  # diversity -------------------------------------------------------------
  dv <- do.call(rbind, lapply(c("Ho", "He", "AR"), function(k) {
    r <- e[e$param %in% paste0(k, "_", pops), ]
    data.frame(pop = factor(plab[sub(paste0("^", k, "_"), "", r$param)], levels = plab), metric = k,
               estimate = r$estimate, lower = r$lower, upper = r$upper)
  }))
  dv$metric <- factor(dv$metric, levels = c("Ho", "He", "AR"),
                      labels = c(fig_label("Ho", "Observed heterozygosity", opts), fig_label("He", "Expected heterozygosity", opts),
                                 fig_label("AR", sprintf("Allelic richness (%d gene copies)", o$comps$rarefaction_genes), opts)))
  figs$diversity <- with_caption(
    ggplot2::ggplot(dv, ggplot2::aes(estimate, pop)) +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.25, colour = eco_pal$ink) +
      ggplot2::facet_wrap(~metric, scales = "free_x", nrow = 1) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::scale_x_continuous(breaks = function(l) { b <- pretty(l, 4); b[b >= l[1] & b <= l[2]] },
                                  expand = ggplot2::expansion(mult = 0.1)) +
      ggplot2::labs(y = lab_pop, x = fig_label("diversity_x", "Estimate (mean over loci)", opts)),
    sprintf("Genetic diversity per population: observed (Ho) and expected (He, Nei's unbiased gene diversity) heterozygosity and allelic richness rarefied to %d gene copies. Points are estimates averaged over loci; bars are %d%% confidence intervals (estimate +/- t(L-1) x SE; SE from a jackknife over loci plus the shared-individual covariance).", o$comps$rarefaction_genes, ci))
  attr(figs$diversity, "fig_size") <- list(width = "double", height = 34 + 6 * length(pops))

  # fst matrix (lower triangle) -----------------------------------------
  r <- length(pops); Fm <- matrix(NA_real_, r, r, dimnames = list(plab, plab))
  pv <- e[match(paste0("Fst_", o$comps$pair_names), e$param), "estimate"]
  Fm[cbind(o$comps$pairs[2, ], o$comps$pairs[1, ])] <- pv
  figs$fst_matrix <- with_caption(
    ps_fst_heat(Fm, fig_label("fst", "Pairwise FST", opts), label = r <= 12),
    sprintf("Pairwise FST (Weir & Cockerham 1984 theta) between sampled populations; lighter = more differentiated. Cell labels give the estimate (leading zero dropped); %d%% confidence intervals for every pair are in the estimates table (Fst_<pop1>_vs_<pop2>).", ci))
  attr(figs$fst_matrix, "fig_size") <- list(width = if (r > 6) "onehalf" else "single", height = if (r > 6) 100 else 72)

  # IBD -------------------------------------------------------------------
  if (!is.null(o$dist) && "ibd_slope" %in% e$param) {
    ok <- o$dist > 0
    y <- pv / (1 - pv)
    lg <- identical(o$ibd_distance, "log")
    x <- if (lg) log(o$dist) else o$dist
    pts <- data.frame(d = o$dist, y = y)[ok, ]
    xs <- seq(min(x[ok]), max(x[ok]), length.out = 60)
    jk <- o$jack[paste0("Fst_", o$comps$pair_names), , drop = FALSE]
    yb <- jk[ok, , drop = FALSE] / (1 - jk[ok, , drop = FALSE])
    xc <- x[ok] - mean(x[ok])
    sl <- colSums(xc * yb) / sum(xc^2); ic <- colMeans(yb) - sl * mean(x[ok])
    pred <- outer(ic, rep(1, length(xs))) + outer(sl, xs)            # L jackknife lines
    Lj <- nrow(pred)
    se_x <- sqrt((Lj - 1) / Lj * colSums(sweep(pred, 2, colMeans(pred))^2))
    infl <- o$se[["ibd_slope"]] / max(o$se_jack[["ibd_slope"]], 1e-12)  # add shared-individual covariance
    tq <- stats::qt(1 - (1 - ci_lvl) / 2, max(L - 1, 1))
    s0 <- sum(xc * y[ok]) / sum(xc^2); i0 <- mean(y[ok]) - s0 * mean(x[ok])
    qq <- rbind(i0 + s0 * xs - tq * infl * se_x, i0 + s0 * xs + tq * infl * se_x)
    band <- data.frame(d = if (lg) exp(xs) else xs, estimate = i0 + s0 * xs, lower = qq[1, ], upper = qq[2, ])
    p <- ggplot2::ggplot(band, ggplot2::aes(d, estimate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.18) +
      ggplot2::geom_line(colour = eco_pal$band, linewidth = 0.8) +
      ggplot2::geom_point(data = pts, ggplot2::aes(d, y), inherit.aes = FALSE, size = 1.4, colour = eco_pal$ink, alpha = 0.8) +
      ggplot2::labs(x = fig_label("distance", if (lg) "Geographic distance (km, log scale)" else "Geographic distance (km)", opts),
                    y = fig_label("ibd_y", "FST / (1 - FST)", opts))
    if (lg) p <- p + ggplot2::scale_x_log10(breaks = function(l) { b <- c(1, 2, 5) * rep(10^(-2:5), each = 3); b[b >= l[1] & b <= l[2]] })
    sr <- e[e$param == "ibd_slope", ]
    figs$ibd <- with_caption(p, sprintf(
      "Isolation by distance: pairwise FST/(1-FST) against %s distance (points = population pairs). Line: least-squares fit (slope %.4f, %d%% CI %.4f to %.4f per %s); band: %d%% confidence band (jackknife over loci plus shared-individual covariance).%s",
      if (lg) "geographic (log axis)" else "geographic", sr$estimate, ci, sr$lower, sr$upper, if (lg) "unit log distance" else "km", ci,
      if (!is.null(o$mantel)) sprintf(" Mantel r = %.2f, P = %.3f.", o$mantel$r, o$mantel$p) else ""))
  }

  # PCA -------------------------------------------------------------------
  if (!is.null(o$pca)) {
    sc <- o$pca$scores
    d <- data.frame(PC1 = sc[, 1], PC2 = sc[, 2], pop = factor(plab[as.character(data$pop)], levels = plab))
    figs$pca <- with_caption(
      ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = pop, shape = pop)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
        ggplot2::geom_point(size = 1.3, alpha = 0.85) +
        ggplot2::scale_colour_manual(values = popcols, breaks = unname(plab), name = lab_pop) +
        ggplot2::scale_shape_manual(values = popshape, breaks = unname(plab), name = lab_pop) +
        ggplot2::labs(x = fig_label("pca_x", sprintf("PC1 (%.1f%% of variance)", o$pca$pct[1]), opts),
                      y = fig_label("pca_y", sprintf("PC2 (%.1f%% of variance)", o$pca$pct[2]), opts)),
      "Principal component analysis of individual genotypes (centred allele frequencies, missing values replaced by the mean). Each point is one individual; colour and symbol show the sampled population; percentages are the share of genetic variance on each axis.")
  }

  # DAPC ------------------------------------------------------------------
  if (!is.null(o$dapc)) {
    dp <- o$dapc
    bic <- dp$bic
    figs$bic_k <- with_caption(
      ggplot2::ggplot(bic, ggplot2::aes(K, BIC)) + ggplot2::geom_line(colour = eco_pal$muted) +
        ggplot2::geom_point(size = 1.4, colour = eco_pal$ink) +
        ggplot2::geom_point(data = bic[bic$K == dp$K, ], size = 3, shape = 21, fill = "#D55E00", colour = eco_pal$ink) +
        ggplot2::scale_x_continuous(breaks = function(l) unique(round(pretty(l)))) +
        ggplot2::labs(x = fig_label("bic_x", "Number of genetic clusters (K)", opts), y = fig_label("bic_y", "BIC", opts)),
      sprintf("BIC of k-means clustering on principal components for K = 1 to %d clusters; the large orange point is the chosen K (%d), from adegenet's '%s' rule. A clear minimum or elbow indicates well-separated clusters; a flat curve indicates weak or continuous structure.",
              max(bic$K), dp$K, o$settings$k_criterion %||% "diffNgroup"))
    if (!is.null(dp$posterior)) {
      post <- dp$posterior
      pop_lab <- factor(plab[as.character(data$pop)], levels = plab)
      ord <- order(data$pop, -apply(post, 1, max), apply(post, 1, which.max))
      cl <- colnames(post)
      cl_lab <- unname(vapply(cl, function(k) fig_label(paste0("cluster_", k), k, opts), ""))
      d <- data.frame(ind = rep(seq_len(nrow(post)), ncol(post)), cluster = factor(rep(cl_lab, each = nrow(post)), levels = cl_lab),
                      p = as.vector(post), pop = rep(pop_lab, ncol(post)))
      d$ind <- factor(d$ind, levels = ord)
      ccol <- stats::setNames(ps_cluster_cols(post, pop_lab, popcols), cl_lab)
      figs$dapc_assignment <- with_caption(
        ggplot2::ggplot(d, ggplot2::aes(ind, p, fill = cluster)) + ggplot2::geom_col(width = 1) +
          ggplot2::facet_grid(~pop, scales = "free_x", space = "free_x") +
          ggplot2::scale_fill_manual(values = ccol, breaks = cl_lab) +
          ggplot2::scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
          ggplot2::scale_x_discrete(breaks = NULL) +
          ggplot2::labs(x = fig_label("dapc_x", "Individuals, grouped by sampled population", opts),
                        y = fig_label("dapc_y", "Membership probability", opts), fill = fig_label("cluster", "Cluster", opts)),
        sprintf("DAPC membership probabilities (one bar per individual, grouped by sampled population) for the K = %d clusters found by k-means, using %d principal component(s) (K - 1; Thia 2023). A cluster made up mainly of one sampled population takes that population's colour (as in the PCA); a cluster shared by several populations is grey.", dp$K, dp$n_pca))
      attr(figs$dapc_assignment, "fig_size") <- list(width = "double", height = 65)
    }
  }

  # map -------------------------------------------------------------------
  if (!is.null(data$coords)) {
    cc <- data$coords[match(pops, data$coords$pop), ]
    cc$He <- e$estimate[match(paste0("He_", pops), e$param)]
    cc$n <- as.integer(table(data$pop)[pops]); cc$lab <- unname(plab)
    # ggrepel point.size is unscaled: mirror scale_size range so labels clear the symbols
    cc$ps <- if (diff(range(cc$n)) > 0) 2 + 3 * (cc$n - min(cc$n)) / diff(range(cc$n)) else 3.5
    ll <- identical(data$coord_type, "lonlat")
    p <- ggplot2::ggplot(cc, ggplot2::aes(x, y)) +
      ggplot2::geom_point(ggplot2::aes(fill = He, size = n), shape = 21, colour = eco_pal$ink, stroke = 0.4) +
      ggplot2::scale_fill_viridis_c(name = fig_label("map_fill", "Expected\nheterozygosity", opts)) +
      ggplot2::scale_size_continuous(name = fig_label("map_size", "Individuals\ngenotyped", opts), range = c(2, 5))
    p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
               ggrepel::geom_text_repel(ggplot2::aes(label = lab, point.size = ps), size = 2.4, colour = eco_pal$ink, seed = 1,
                                        box.padding = 0.25, point.padding = 0.1, min.segment.length = 0.2,
                                        segment.size = 0.25, segment.colour = eco_pal$muted)
             else ggplot2::geom_text(ggplot2::aes(label = lab), size = 2.4, vjust = -1.4)
    p <- p + ggplot2::labs(x = fig_label("map_x", if (ll) "Longitude (°)" else "Easting (km)", opts),
                           y = fig_label("map_y", if (ll) "Latitude (°)" else "Northing (km)", opts)) +
      ggplot2::scale_x_continuous(n.breaks = 4, expand = ggplot2::expansion(mult = 0.15)) +
      ggplot2::scale_y_continuous(n.breaks = 5, expand = ggplot2::expansion(mult = 0.1))
    p <- p + if (ll) ggplot2::coord_quickmap() else ggplot2::coord_equal()
    figs$sample_map <- with_caption(p, sprintf(
      "Sampled populations (%s): point fill = expected heterozygosity (He, Nei's unbiased gene diversity), point size = number of individuals genotyped. He intervals are in the diversity figure.",
      if (ll) "longitude/latitude in decimal degrees" else "coordinates in km"))
    attr(figs$sample_map, "fig_size") <- list(width = "single", height = 95)
  }
  figs
}

# Pairwise FST heatmap, lower triangle only (empty cells dropped, not grey).
# fig_matrix() draws the full square with NA cells; the triangle avoids
# repeating every value twice. Label ink switches by position on the viridis
# scale (white on dark, black on light), not by the median.
ps_fst_heat <- function(Fm, fill_lab, label = TRUE) {
  d <- as.data.frame(as.table(Fm)); names(d) <- c("row", "col", "value")
  d <- d[!is.na(d$value), ]
  lv <- rownames(Fm)
  d$row <- factor(d$row, levels = rev(lv)); d$col <- factor(d$col, levels = lv)
  lim <- c(min(0, d$value), max(d$value))
  p <- ggplot2::ggplot(d, ggplot2::aes(col, row, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) + ggplot2::coord_equal() +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
    scale_fill_eco_seq(limits = lim) +
    ggplot2::labs(x = NULL, y = NULL, fill = fill_lab)
  if (label) {
    d$light <- (d$value - lim[1]) / max(diff(lim), 1e-12) > 0.55
    p <- p + ggplot2::geom_text(data = d, ggplot2::aes(label = sub("^(-?)0\\.", "\\1.", sprintf("%.2f", value)), colour = light),
                                size = 2.1, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none")
  }
  p
}

# Slope interval (resampling loci and individuals) and Mantel test (permutation of populations) disagree.
ps_ibd_mixed <- function(e, o, mt) {
  pops <- o$comps$pops; pr <- o$comps$pairs
  f <- e$estimate[match(paste0("Fst_", o$comps$pair_names), e$param)]
  mean_f <- vapply(seq_along(pops), function(i) mean(f[pr[1, ] == i | pr[2, ] == i]), 1)
  out <- pops[which.max(mean_f)]
  sprintf(paste("The evidence is mixed. The slope interval, which reflects uncertainty from the sample of loci and individuals for these particular populations, excludes 0, but the Mantel test, which asks whether a pattern this strong could arise from an arbitrary arrangement of these populations in space, does not support it (P = %.2f).",
                "%s is strongly differentiated from all others (mean pairwise FST %.2f vs %.2f overall), regardless of distance; one small or isolated population like this inflates or masks the distance relationship.",
                "Do not conclude that differentiation increases with distance; re-run without %s (or with more populations) to check."),
          mt$p, out, max(mean_f), mean(f), out)
}
