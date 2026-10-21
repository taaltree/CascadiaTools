# Interpretation, figures and methods text for popgen-ne.

pg_references <- c(
  "Waples RS (2006) A bias correction for estimates of effective population size based on linkage disequilibrium at unlinked gene loci. Conservation Genetics 7:167-184.",
  "Waples RS, Do C (2008) LDNE: a program for estimating effective population size from data on linkage disequilibrium. Molecular Ecology Resources 8:753-756.",
  "Waples RS, Do C (2010) Linkage disequilibrium estimates of contemporary Ne using highly variable genetic markers: a largely untapped resource for applied conservation and evolution. Evolutionary Applications 3:244-262.",
  "Peel D, Waples RS, Macbeth GM, Do C, Ovenden JR (2013) Accounting for missing data in the estimation of contemporary genetic effective population size (Ne). Molecular Ecology Resources 13:243-253.",
  "Do C, Waples RS, Peel D, Macbeth GM, Tillett BJ, Ovenden JR (2014) NeEstimator v2: re-implementation of software for the estimation of contemporary effective population size (Ne) from genetic data. Molecular Ecology Resources 14:209-214.",
  "Jones AT, Ovenden JR, Wang Y-G (2016) Improved confidence intervals for the linkage disequilibrium method for estimating effective population size. Heredity 117:217-223.",
  "Nei M, Tajima F (1981) Genetic drift and estimation of effective population size. Genetics 98:625-640.",
  "Pollak E (1983) A new method for estimating the effective population size from allele frequency changes. Genetics 104:531-548.",
  "Waples RS (1989) A generalized approach for estimating effective population size from temporal changes in allele frequency. Genetics 121:379-391.",
  "Jorde PE, Ryman N (2007) Unbiased estimator for genetic drift and effective population size. Genetics 177:927-935.",
  "Queller DC, Goodnight KF (1989) Estimating relatedness using genetic markers. Evolution 43:258-275.",
  "Lynch M, Ritland K (1999) Estimation of pairwise relatedness with molecular markers. Genetics 152:1753-1766.",
  "Wang J (2014) Marker-based estimates of relatedness and inbreeding coefficients: an assessment of current methods. Journal of Evolutionary Biology 27:518-530.",
  "Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen W-M (2010) Robust relationship inference in genome-wide association studies. Bioinformatics 26:2867-2873.",
  "Cornuet JM, Luikart G (1996) Description and power analysis of two tests for detecting recent population bottlenecks from allele frequency data. Genetics 144:2001-2014.",
  "Piry S, Luikart G, Cornuet J-M (1999) BOTTLENECK: a computer program for detecting recent reductions in the effective size using allele frequency data. Journal of Heredity 90:502-503.",
  "Franklin IR (1980) Evolutionary change in small populations. In: Soule ME, Wilcox BA (eds) Conservation Biology: An Evolutionary-Ecological Perspective. Sinauer, pp 135-150.",
  "Frankham R, Bradshaw CJA, Brook BW (2014) Genetics in conservation management: revised recommendations for the 50/500 rules, Red List criteria and population viability analyses. Biological Conservation 170:56-63.",
  "Franklin IR, Allendorf FW, Jamieson IG (2014) The 50/500 rule is still valid - reply to Frankham et al. Biological Conservation (doi:10.1016/j.biocon.2014.05.004).",
  "Frankham R (1995) Effective population size/adult population size ratios in wildlife: a review. Genetical Research 66:95-107.",
  "Waples RS, Anderson EC (2017) Purging putative siblings from population genetic data sets: a cautionary view. Molecular Ecology 26:1211-1224.",
  "Waples RS, England PR (2011) Estimating contemporary effective population size on the basis of linkage disequilibrium in the face of migration. Genetics 189:633-644.",
  "Waples RS, Antao T, Luikart G (2014) Effects of overlapping generations on linkage disequilibrium estimates of effective population size. Genetics 197:769-780."
)

fmt_ne <- function(x) ifelse(!is.finite(x) | x <= 0, "infinity", formatC(x, format = "f", digits = if (abs(x) < 100) 1 else 0, big.mark = ","))
fmt_ci <- function(lo, hi) sprintf("%s to %s", fmt_ne(lo), fmt_ne(hi))

ne_category <- function(ne, lo, hi) {
  if (!is.finite(hi) || hi <= 0) hi <- Inf
  if (hi < 50) return("The whole 95% interval lies below 50, the classic short-term threshold (Franklin 1980) and below the revised 100 (Frankham et al. 2014): inbreeding is accumulating fast enough that inbreeding depression is a real short-term risk. Genetic rescue or increased connectivity deserve consideration.")
  if (hi < 100) return("The whole 95% interval lies below 100, the revised short-term threshold of Frankham et al. (2014) (the older 50/500 rule used 50): the population is small enough that inbreeding depression is a concern over the next few generations.")
  if (lo >= 1000) return("The 95% interval lies entirely above 1,000, the revised long-term threshold for retaining evolutionary potential (Frankham et al. 2014; 500 in the older rule).")
  if (lo >= 500 && hi >= 1000) return("The lower 95% limit exceeds 500 (the classic long-term threshold, Franklin 1980) but the interval reaches below 1,000 (the revised threshold of Frankham et al. 2014): short-term inbreeding is not a concern; long-term adaptive potential is probably but not certainly adequate.")
  if (lo >= 100) return("The lower 95% limit exceeds 100, so short-term inbreeding depression is unlikely to be a concern (Frankham et al. 2014), but Ne may be below the 500-1,000 needed to retain evolutionary potential in the long term.")
  if (lo >= 50) return("The interval lies above 50 but reaches below 100: the population passes the classic short-term threshold (50) but not clearly the revised one (100; Frankham et al. 2014).")
  "The 95% interval spans the 50/100 short-term thresholds, so these data cannot say whether the population is above or below them; more individuals or loci, or a second sample, would narrow it."
}

pg_interpret <- function(fit, data) {
  ob <- fit$object; st <- ob$settings
  est <- fit$estimates
  out <- character()
  for (s in names(ob$ld)) {
    r <- ob$ld[[s]]
    if (is.na(r$inv)) next
    h <- r$inv / 2
    h10 <- 1 - (1 - max(h, 0))^10
    txt <- sprintf("Sample %s (%d individuals, %d loci after excluding alleles rarer than %.3f): the linkage-disequilibrium estimate of Ne is %s (95%% CI %s; %s).",
                   s, r$n_ind, r$n_loci, r$pcrit, fmt_ne(r$ne), fmt_ci(r$ne_lo, r$ne_hi), r$ci_type)
    if (is.finite(r$ne) && r$ne > 0)
      txt <- paste(txt, sprintf("An Ne of this size loses about %.2f%% of its heterozygosity per generation (1/(2Ne)), about %.1f%% over 10 generations, unless immigration or mutation replaces it.",
                                100 * h, 100 * h10))
    else txt <- paste(txt, "The point estimate is infinite: the LD signal of drift is no larger than expected from sampling alone. This happens when the true Ne is large relative to the number of individuals sampled; read the lower limit as the useful information ('Ne is probably at least", fmt_ne(r$ne_lo), "').")
    if (!is.finite(r$ne_hi) && is.finite(r$ne)) txt <- paste(txt, sprintf("The upper limit is infinite: the data rule out an Ne below %s but cannot exclude a very large population.", fmt_ne(r$ne_lo)))
    out <- c(out, txt, ne_category(r$ne, r$ne_lo, r$ne_hi))
    if (!is.null(r$by_pcrit) && nrow(r$by_pcrit) > 1) {
      fin <- r$by_pcrit$ne[is.finite(r$by_pcrit$ne)]
      if (length(fin) > 1)
        out <- c(out, sprintf("Sensitivity to the minimum allele frequency (Pcrit %s): point estimates range from %s to %s. Waples & Do (2010) recommend Pcrit = 0.02 for samples above 25 individuals; very low Pcrit admits rare alleles that bias LD estimates upward.",
                              paste(r$by_pcrit$pcrit, collapse = ", "), fmt_ne(min(fin)), fmt_ne(max(fin))))
    }
  }
  for (pp in names(ob$temporal)) {
    tm <- ob$temporal[[pp]]
    suf <- if (length(ob$temporal) == 1) "" else paste0("_", pg_label(pp))
    e <- est[est$param == paste0("Ne_temporal", suf), ]
    out <- c(out, sprintf("Temporal method (population %s, samples at times %s and %s, %.1f generations apart, sampling plan %s%s): Jorde & Ryman's Fs' = %.4f gives Ne = %s (95%% CI %s; %s). Nei-Tajima and Pollak estimators give %s and %s. The temporal estimate is a harmonic mean of Ne over the generations between the samples; the LD estimate mostly reflects the generation(s) just before the later sample.",
                          pp, tm$times[1], tm$times[2], tm$t, tm$plan,
                          if (tm$plan == "I") sprintf(", census N = %s", format(tm$census_N)) else "",
                          tm$JR$Fp, fmt_ne(e$estimate), fmt_ci(e$lower, e$upper), tm$JR$ci_type,
                          fmt_ne(est$estimate[est$param == paste0("Ne_temporal_NT", suf)]),
                          fmt_ne(est$estimate[est$param == paste0("Ne_temporal_Pollak", suf)])))
    if (tm$t < 2) out <- c(out, "With fewer than 2 generations between samples, the drift signal is small compared with sampling noise and the temporal estimate is imprecise.")
  }
  for (s in names(ob$rel)) {
    rr <- ob$rel[[s]]; p <- rr$pairs
    th <- st$kin_thresholds
    n1 <- sum(p$r >= th[1], na.rm = TRUE); n2 <- sum(p$r >= th[2] & p$r < th[1], na.rm = TRUE)
    txt <- sprintf("Relatedness in sample %s (%d pairs, %d loci, %s estimator): %d pairs have r >= %.3f (likely first-degree relatives: full sibs or parent-offspring) and %d have %.3f <= r < %.3f (likely second-degree, e.g. half sibs). Mean pairwise r is %.3f; moment estimators measure relatedness relative to this sample's allele frequencies, so unrelated pairs average slightly below zero.",
                   s, nrow(p), rr$n_loci, if (st$rel_primary == "LR") "Lynch-Ritland" else "Queller-Goodnight", n1, th[1], n2, th[2], th[1], mean(p$r, na.rm = TRUE))
    if (!is.null(rr$error)) {
      e <- rr$error; eu <- e[e$class == "unrelated", ]; ef <- e[e$class == "fullsib", ]; eh <- e[e$class == "halfsib", ]
      nu <- nrow(p)
      txt <- paste(txt, sprintf("Simulated pairs drawn from these allele frequencies show what the markers can resolve: unrelated pairs score r = %.3f +/- %.3f (SD), so about %.1f%% of truly unrelated pairs would be called second-degree or closer (~%.0f false calls among %d pairs); %.0f%% of full sibs and %.0f%% of half sibs fall in their own class. These simulated rates assume unrelated founders; in a small population real half sibs score lower (background kinship, sample-based allele frequencies), so in validation runs more true half sibs fell below the second-degree threshold than this simulation predicts (10%% vs 2%% with 1,000 SNPs; 40%% vs 34%% with 15 microsatellites).",
                                eu$mean_r, eu$sd_r, 100 * (eu$p_first + eu$p_second), nu * (eu$p_first + eu$p_second), nu, 100 * ef$p_first, 100 * eh$p_second))
    }
    out <- c(out, txt)
    if (!is.null(rr$known)) for (cl in names(rr$known)) {
      k <- rr$known[[cl]]; k1 <- k[k$estimator == st$rel_primary, ]
      expv <- c(fullsib = 0.5, parent_offspring = 0.5, parentoffspring = 0.5, halfsib = 0.25, unrelated = 0)[cl]
      out <- c(out, sprintf("Known %s pairs (n = %d) have mean r = %.3f (95%% CI %.3f to %.3f)%s.", cl, k1$n_pairs, k1$estimate, k1$lower, k1$upper,
                            if (!is.na(expv)) sprintf(", against a pedigree expectation of %.2f", expv) else ""))
    }
  }
  for (s in names(ob$het)) {
    h <- ob$het[[s]]$table
    out <- c(out, sprintf("Sample %s: observed heterozygosity %.3f, expected %.3f, F_IS = %.3f (95%% CI %.3f to %.3f)%s.", s,
                          h$estimate[1], h$estimate[2], h$estimate[3], h$lower[3], h$upper[3],
                          if (h$lower[3] > 0) "; heterozygote deficit (possible inbreeding, null alleles or a mixture of subpopulations)"
                          else if (h$upper[3] < 0) "; heterozygote excess (expected when few parents produced the sample)" else "; consistent with random mating"))
    b <- ob$het[[s]]$bottleneck
    if (!is.null(b)) out <- c(out, sprintf("Bottleneck test (heterozygosity excess, infinite-alleles model; Cornuet & Luikart 1996): %d of %d loci show He above the mutation-drift equilibrium expectation; Monte Carlo P = %.3g for an excess of heterozygosity (mean standardised difference over loci compared with its equilibrium distribution; the BOTTLENECK-style Wilcoxon test, P = %.3g, is shown for reference only because it is anti-conservative here). A small P is consistent with a recent bottleneck, but also with any history that leaves the population out of mutation-drift equilibrium (e.g. recent admixture or founding); this test assumes the population was at mutation-drift equilibrium before any bottleneck and is sensitive to the mutation model.",
                                           b$n_excess, b$n_loci, b$p_value, b$p_wilcoxon))
  }
  out <- c(out, "The thresholds are rules of thumb, not cliffs: Franklin (1980) proposed Ne = 50 to keep inbreeding at about 1% per generation in the short term and 500 to retain adaptive potential; Frankham et al. (2014) argued these should be raised to 100 and 1,000, while Franklin, Allendorf & Jamieson (2014) defended the original values. Whichever is used, compare the whole confidence interval, not only the point estimate, with the threshold, and treat a single estimate as one line of evidence alongside census trends and connectivity.")
  out <- c(out, "Ne is not census size: Ne/N averages roughly 0.1-0.2 in wild populations (Frankham 1995), so an Ne of 100 typically corresponds to several hundred adults. In species with overlapping generations a single-cohort LD estimate approaches the number of breeders in one season rather than per-generation Ne (Waples et al. 2014); immigration and sampling of several subpopulations bias LD estimates (Waples & England 2011), and over-representing a few families biases them downward (Waples & Anderson 2017). In ecotools' independent validation with separate sexes (no selfing), LD estimates ran about 5-25% above the true Ne (harmonic mean over replicates; the larger values with 20 microsatellites, a skewed sex ratio or a sample that is a large share of the population), whereas the temporal estimate was within 6%; treat a single LD estimate as slightly optimistic.")
  out
}

pg_methods <- function(fit, data, settings) {
  ob <- fit$object; st <- ob$settings
  pc <- vapply(ob$ld, function(r) r$pcrit, 0)
  sprintf(paste(
    "Effective population size was estimated with the linkage-disequilibrium method (LDNe; Waples 2006; Waples & Do 2008), implemented in R following the NeEstimator v2 formulas (Do et al. 2014): Burrows' composite r^2 between all pairs of loci (S/(S-1) adjusted), averaged over allele pairs and weighted by the number of independent allelic comparisons (and by sample size squared when data were missing; Peel et al. 2013), Waples' (2006) sample-size bias correction, a %s mating model, and exclusion of alleles with frequency below Pcrit = %s (Waples & Do 2010). 95%% confidence intervals were obtained by jackknifing over individuals with the 0.84 variance correction of Jones et al. (2016).",
    "Where a population was sampled twice, temporal Ne was estimated from allele-frequency change with Jorde & Ryman's (2007) Fs' (95%% CI: chi-square interval on the uncorrected Fs with degrees of freedom 2/CV^2 estimated by a jackknife over loci, the analogue of the LD jackknife CI, mapped to Fs'; the jackknife variance is inflated by a design effect 1 + 0.5 (D - D/L) R^2' for the correlation of allele-frequency change among unlinked loci in drift linkage disequilibrium, where D is the number of independent alleles, L the number of loci and R^2' the drift LD from the LD method; calibrated by simulation in the ecotools validation), with Nei & Tajima (1981) Fc and Pollak (1983) Fk (chi-square CIs) for comparison, using sampling plan %s (Waples 1989).",
    "Pairwise relatedness was estimated with the Lynch & Ritland (1999) and Queller & Goodnight (1989) moment estimators using within-sample allele frequencies; pairs were classed as first-degree (r >= %.3f) or second-degree (r >= %.3f) relatives following the KING thresholds (Manichaikul et al. 2010), and misclassification rates at these markers were estimated by simulating unrelated, half-sib, full-sib and parent-offspring pairs from the observed allele frequencies.",
    "Observed and expected heterozygosity and F_IS were averaged over loci with jackknife CIs. The method was validated by simulation (ecotools popgen-ne VALIDATION.md)."),
    st$mating, paste(unique(signif(pc, 3)), collapse = "/"),
    if (length(ob$temporal)) paste(unique(vapply(ob$temporal, function(x) x$plan, "")), collapse = "/") else "II",
    st$kin_thresholds[1], st$kin_thresholds[2])
}

# --- figures -----------------------------------------------------------------
# Publication suite (docs/FIGURES.md): ne_estimates, ne_by_maf,
# relatedness_distribution, kin_resolution, kin_network, design_precision.
# Every label goes through fig_label(); colours come from opts$palette in a
# fixed entity order (methods, samples, relationship classes), so a class keeps
# its colour in every figure and greyscale works (shape/linetype as a second
# encoding). Intervals are the ones computed at fit time (settings: ci);
# frequentist, so there is no inner band.

# Fixed colours for an ordered set of entities. With Okabe-Ito the yellow
# (4th) is skipped for lines/points on white; other palettes are used in order.
pg_entity_cols <- function(levels, opts) {
  pal <- opts$palette %||% "okabe_ito"
  full <- if (length(pal) > 1) pal else eco_pal[[pal]]
  if (is.null(full)) full <- eco_pal$okabe_ito
  if (identical(pal, "okabe_ito")) full <- full[-4]
  if (identical(pal, "grey")) full <- c("#1a1a1a", "#595959", "#8c8c8c", "#b0b0b0")   # lightest core grey vanishes on white
  if (length(levels) > length(full)) full <- rep(full, length.out = length(levels))
  stats::setNames(full[seq_along(levels)], levels)
}
# Size hint that yields to a user's global `figures: width:` (core lets a module
# hint beat the global width; per_figure still beats both).
pg_size <- function(width, height, opts) {
  w <- if (!is.null(opts$width) && !identical(opts$width, "single")) opts$width else width
  list(width = w, height = height)
}
# Effective width in mm of one figure (per_figure > global non-default > hint).
pg_width_mm <- function(name, hint, opts) {
  w <- opts$per_figure[[name]]$width %||% pg_size(hint, NA, opts)$width
  wm <- tryCatch(fig_width_mm(w, opts), error = function(e) NULL)
  if (is.null(wm)) 120 else wm
}
# Relationship classes: canonical keys (as io.R normalises them) and labels.
pg_rel_keys <- c("unrelated", "halfsib", "fullsib", "parent_offspring")
pg_rel_canon <- function(x) { x <- sub("^parentoffspring$", "parent_offspring", x); x }
pg_rel_labels <- function(opts) c(
  unrelated = fig_label("unrelated", "Unrelated", opts), halfsib = fig_label("halfsib", "Half-sib", opts),
  fullsib = fig_label("fullsib", "Full-sib", opts), parent_offspring = fig_label("parent_offspring", "Parent-offspring", opts))

pg_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  ob <- fit$object; st <- ob$settings; est <- fit$estimates
  ci <- if (is.null(st$ci)) 0.95 else st$ci          # level the intervals were computed at
  if (!is.null(opts$ci) && abs(opts$ci - ci) > 1e-9)
    message("popgen-ne figures: intervals are ", round(100 * ci), "% (computed at fit time); set settings: ci: ",
            opts$ci, " and refit to draw ", round(100 * opts$ci), "% intervals.")
  pct <- round(100 * ci)
  inf_lab <- fig_label("infinite", "infinite", opts)
  rel_lab <- pg_rel_labels(opts)
  rel_col <- pg_entity_cols(rel_lab[pg_rel_keys], opts)
  rel_lty <- stats::setNames(c("dotted", "dashed", "solid", "dotdash"), rel_lab[pg_rel_keys])
  grey <- identical(opts$palette, "grey") || isTRUE(opts$greyscale)
  thr_col <- if (grey) eco_pal$ink else eco_pal$truth        # kin thresholds / study design marker
  thr_word <- if (grey) "black" else "red"
  figs <- list()

  # ne_estimates: all Ne rows (LD per sample, temporal per pop and method)
  ne <- est[grepl("^Ne_", est$param), ]
  if (nrow(ne)) {
    cap <- 1e5
    mkey <- ifelse(grepl("^Ne_LD", ne$param), "LD", ifelse(grepl("_NT", ne$param), "NT", ifelse(grepl("_Pollak", ne$param), "Pollak", "JR")))
    mlab <- c(LD = fig_label("method_LD", "LD (single sample)", opts), JR = fig_label("method_JR", "Temporal: Jorde-Ryman", opts),
              NT = fig_label("method_NT", "Temporal: Nei-Tajima", opts), Pollak = fig_label("method_Pollak", "Temporal: Pollak", opts))
    mshort <- c(LD = "LD", JR = "Jorde-Ryman", NT = "Nei-Tajima", Pollak = "Pollak")
    d <- data.frame(param = ne$param, method = factor(mlab[mkey], levels = mlab),
                    estimate = pmin(ifelse(is.finite(ne$estimate) & ne$estimate > 0, ne$estimate, cap), cap),
                    lower = pmin(ifelse(is.finite(ne$lower) & ne$lower > 0, ne$lower, 1), cap),
                    upper = pmin(ifelse(is.finite(ne$upper) & ne$upper > 0, ne$upper, cap), cap),
                    open = !is.finite(ne$upper) | ne$upper <= 0 | ne$upper > cap)
    tpops <- names(ob$temporal)
    d$label <- vapply(d$param, function(p) {
      if (grepl("^Ne_LD_", p)) { s <- sub("^Ne_LD_", "", p); return(fig_label(s, gsub("_", " ", s), opts)) }
      pp <- sub("^Ne_temporal(_NT|_Pollak)?_?", "", p)
      tm <- if (nzchar(pp)) ob$temporal[[which(pg_label(tpops) == pp)[1]]] else ob$temporal[[1]]
      paste0(if (length(tpops) > 1) paste0(fig_label(tm$pop, tm$pop, opts), " ") else "", tm$times[1], "-", tm$times[2])
    }, "")
    d$label <- paste0(d$label, ": ", vapply(mkey, function(k) fig_label(paste0("short_", k), mshort[[k]], opts), ""))
    d$label <- vapply(d$label, function(l) fig_label(l, l, opts), "")
    d$label <- factor(d$label, levels = rev(unique(d$label)))
    shp <- stats::setNames(c(16, 17, 15, 18), mlab)
    narrow <- pg_width_mm("ne_estimates", "onehalf", opts) < 100   # long row labels squeeze a single-column panel
    p <- ggplot2::ggplot(d, ggplot2::aes(estimate, label, colour = method, shape = method)) +
      ggplot2::geom_vline(xintercept = c(50, 500), linetype = 2, colour = eco_pal$null, linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = c(100, 1000), linetype = 3, colour = eco_pal$null, linewidth = 0.35) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.6) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_point(data = d[d$open, ], ggplot2::aes(x = upper), shape = 62, size = 3.5, show.legend = FALSE) +
      ggplot2::scale_x_log10(breaks = c(10, 50, 100, 500, 1000, 1e4, 1e5), labels = c("10", "50", "100", "500", "1,000", "10,000", inf_lab),
                             expand = ggplot2::expansion(mult = c(0.04, 0.08))) +
      ggplot2::scale_y_discrete(labels = function(x) if (narrow) sub(": ", ":\n", x, fixed = TRUE) else x) +
      ggplot2::scale_colour_manual(values = pg_entity_cols(mlab, opts), drop = TRUE) +
      ggplot2::scale_shape_manual(values = shp, drop = TRUE) +
      ggplot2::labs(x = fig_label("Ne", "Effective population size Ne (log scale)", opts), y = NULL, colour = NULL, shape = NULL) +
      ggplot2::guides(colour = "none", shape = "none")
    attr(p, "fig_size") <- pg_size("onehalf", 25 + (if (narrow) 11 else 8) * length(unique(d$label)), opts)
    figs$ne_estimates <- with_caption(p, sprintf("Effective population size by sample and method, with %d%% confidence intervals (LD: jackknife over individuals; Jorde-Ryman: jackknife over loci; Nei-Tajima and Pollak: chi-square). Colour and symbol show the method. Dashed lines: 50 and 500 (Franklin 1980); dotted lines: the revised 100 and 1,000 (Frankham et al. 2014). An arrow (>) marks an infinite upper limit; infinite point estimates are drawn at the right edge.", pct))
  }

  # ne_by_maf: LD estimate against the allele-frequency cutoff, per sample
  bp <- do.call(rbind, lapply(names(ob$ld), function(s) if (!is.null(ob$ld[[s]]$by_pcrit)) data.frame(sample = s, ob$ld[[s]]$by_pcrit)))
  if (!is.null(bp) && nrow(bp)) {
    cap <- 1e5
    bp$est <- pmin(ifelse(is.finite(bp$ne) & bp$ne > 0, bp$ne, cap), cap)
    bp$lo <- pmin(ifelse(is.finite(bp$ne_lo) & bp$ne_lo > 0, bp$ne_lo, 1), cap)
    bp$hi <- pmin(ifelse(is.finite(bp$ne_hi) & bp$ne_hi > 0, bp$ne_hi, cap), cap)
    slev <- vapply(names(ob$ld), function(s) fig_label(s, gsub("_", " ", s), opts), "")
    bp$sample <- factor(slev[bp$sample], levels = unique(slev))
    pd <- ggplot2::position_dodge(width = 0.35)
    bp$pc <- factor(format(bp$pcrit), levels = format(sort(unique(bp$pcrit))))
    p <- ggplot2::ggplot(bp, ggplot2::aes(pc, est, colour = sample, shape = sample, linetype = sample, group = sample)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), position = pd, linewidth = 0.5, linetype = 1) +
      ggplot2::geom_point(position = pd, size = 1.8) + ggplot2::geom_line(position = pd, linewidth = 0.4) +
      ggplot2::scale_y_log10(breaks = c(10, 30, 100, 300, 1000, 3000, 1e4, 1e5), labels = c("10", "30", "100", "300", "1,000", "3,000", "10,000", inf_lab)) +
      ggplot2::scale_colour_manual(values = pg_entity_cols(levels(bp$sample), opts)) +
      ggplot2::labs(x = fig_label("pcrit", "Minimum allele frequency used, Pcrit", opts),
                    y = fig_label("Ne_LD", "LD estimate of Ne (log scale)", opts), colour = NULL, shape = NULL, linetype = NULL)
    if (length(levels(bp$sample)) > 1) p <- p + ggplot2::guides(colour = ggplot2::guide_legend(position = "bottom"), shape = ggplot2::guide_legend(position = "bottom"), linetype = ggplot2::guide_legend(position = "bottom"))
    # core honours a height hint (not an aspect hint) for module figures
    attr(p, "fig_size") <- pg_size("onehalf", 85, opts)
    figs$ne_by_maf <- with_caption(p, sprintf("Sensitivity of the LD estimate of Ne to the minimum allele frequency Pcrit, with %d%% jackknife confidence intervals (points: estimates; one line per sample). Rare alleles inflate LD estimates when Pcrit is too low (Waples & Do 2010). Values at the top ('%s') have an infinite estimate or limit.", pct, inf_lab))
  }

  # relatedness_distribution: histogram of all pairs + known pairs as a rug
  rel <- do.call(rbind, lapply(names(ob$rel), function(s) { p <- ob$rel[[s]]$pairs; data.frame(sample = s, r = p$r, known = if (is.null(p$known)) NA else p$known) }))
  rel_xlab <- fig_label("relatedness", sprintf("Pairwise relatedness r (%s)", if (identical(st$rel_primary, "QG")) "Queller-Goodnight" else "Lynch-Ritland"), opts)
  if (!is.null(rel) && nrow(rel)) {
    rel$sample <- vapply(rel$sample, function(s) fig_label(s, gsub("_", " ", s), opts), "")
    refs <- data.frame(x = c(0, 0.25, 0.5))
    p <- ggplot2::ggplot(rel[!is.na(rel$r) & is.na(rel$known), ], ggplot2::aes(r)) +
      ggplot2::geom_histogram(bins = 60, fill = eco_pal$grey[3], colour = "white", linewidth = 0.2) +
      ggplot2::geom_vline(data = refs, ggplot2::aes(xintercept = x), linetype = 2, colour = eco_pal$null, linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = st$kin_thresholds, colour = thr_col, linewidth = 0.4) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.05))) +
      ggplot2::labs(x = rel_xlab, y = fig_label("pairs", "Number of pairs", opts))
    kn <- rel[!is.na(rel$known) & !is.na(rel$r), ]
    if (nrow(kn)) {
      kn$known <- pg_rel_canon(kn$known)
      lab <- ifelse(kn$known %in% names(rel_lab), rel_lab[kn$known], vapply(kn$known, function(k) fig_label(k, k, opts), ""))
      lv <- c(rel_lab[pg_rel_keys], setdiff(unique(lab), rel_lab))
      kn$group <- factor(lab, levels = lv[lv %in% lab])
      cols <- c(rel_col, pg_entity_cols(setdiff(levels(kn$group), names(rel_col)), opts))
      p <- p + ggplot2::geom_rug(data = kn, ggplot2::aes(x = r, colour = group, linetype = group), sides = "b", length = ggplot2::unit(0.06, "npc"), linewidth = 0.6) +
        ggplot2::scale_colour_manual(values = cols, drop = TRUE) +
        ggplot2::scale_linetype_manual(values = c(rel_lty, stats::setNames(rep("solid", length(setdiff(levels(kn$group), names(rel_lty)))), setdiff(levels(kn$group), names(rel_lty)))), drop = TRUE) +
        ggplot2::labs(colour = fig_label("known_relationship", "Known relationship", opts), linetype = fig_label("known_relationship", "Known relationship", opts)) +
        ggplot2::guides(colour = ggplot2::guide_legend(position = "top"), linetype = ggplot2::guide_legend(position = "top"))
    }
    if (length(unique(rel$sample)) > 1) { p <- p + ggplot2::facet_wrap(~sample, scales = "free_y"); attr(p, "fig_size") <- pg_size("double", 85, opts) }
    figs$relatedness_distribution <- with_caption(p, sprintf("Distribution of pairwise relatedness among sampled pairs of unknown relationship (bars)%s. Grey dashed lines: pedigree expectations for unrelated (0), half-sib (0.25) and full-sib/parent-offspring (0.5) pairs; %s lines: thresholds used to call second- (r >= %.3f) and first-degree (r >= %.3f) relatives.%s",
      if (length(unique(rel$sample)) > 1) ", one panel per sample" else "", if (grey) "Solid black" else "Red", st$kin_thresholds[2], st$kin_thresholds[1],
      if (nrow(kn)) sprintf(" Ticks below the axis: pairs of known relationship (%s).", if (length(levels(kn$group)) > 1) "colour and line type show the class" else levels(kn$group)) else ""))

    # kin_resolution: simulated pairs of each class at these markers
    simd <- do.call(rbind, lapply(names(ob$rel), function(s) if (!is.null(ob$rel[[s]]$sim)) data.frame(sample = s, ob$rel[[s]]$sim)))
    if (!is.null(simd) && nrow(simd)) {
      simd$sample <- vapply(simd$sample, function(s) fig_label(s, gsub("_", " ", s), opts), "")
      simd$class <- factor(rel_lab[pg_rel_canon(simd$class)], levels = rel_lab[pg_rel_keys])
      simd <- simd[!is.na(simd$class) & is.finite(simd$r), ]
      p2 <- ggplot2::ggplot(simd, ggplot2::aes(r, fill = class, colour = class, linetype = class)) +
        ggplot2::geom_density(alpha = 0.35, linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = st$kin_thresholds, colour = thr_col, linewidth = 0.4) +
        ggplot2::scale_fill_manual(values = rel_col, drop = TRUE) + ggplot2::scale_colour_manual(values = rel_col, drop = TRUE) +
        ggplot2::scale_linetype_manual(values = rel_lty, drop = TRUE) +
        ggplot2::scale_x_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
        ggplot2::labs(x = fig_label("relatedness_sim", sprintf("Estimated r of simulated pairs (%s)", if (identical(st$rel_primary, "QG")) "Queller-Goodnight" else "Lynch-Ritland"), opts),
                      y = fig_label("density", "Probability density", opts), fill = NULL, colour = NULL, linetype = NULL) +
        ggplot2::guides(fill = ggplot2::guide_legend(position = "top"), colour = ggplot2::guide_legend(position = "top"), linetype = ggplot2::guide_legend(position = "top"))
      if (length(unique(simd$sample)) > 1) { p2 <- p2 + ggplot2::facet_wrap(~sample); attr(p2, "fig_size") <- pg_size("double", 85, opts) }
      figs$kin_resolution <- with_caption(p2, sprintf("How well these markers separate relationship classes: densities of relatedness estimates for pairs simulated from the sample's allele frequencies (unlinked loci, Hardy-Weinberg), by true class. Overlap between classes is the misclassification rate; %s lines are the calling thresholds (%.3f and %.3f).", thr_word, st$kin_thresholds[2], st$kin_thresholds[1]))
    }

    # kin_network: pairs above the second-degree threshold (first-degree only if too many)
    kp <- do.call(rbind, lapply(names(ob$rel), function(s) { p <- ob$rel[[s]]$pairs; p <- p[!is.na(p$r) & p$r >= st$kin_thresholds[2], ]; if (nrow(p)) data.frame(sample = s, p) }))
    only_first <- FALSE
    if (!is.null(kp) && nrow(kp) > 120) { kp <- kp[kp$r >= st$kin_thresholds[1], ]; only_first <- TRUE }
    if (!is.null(kp) && nrow(kp)) {
      top <- nrow(kp) > 300
      if (top) kp <- kp[order(-kp$r)[1:300], ]
      kp <- kp[order(kp$id1), ]
      nodes <- unique(as.vector(rbind(kp$id1, kp$id2))); nn <- length(nodes)   # partners placed next to each other
      ang <- pi / 2 - 2 * pi * (seq_len(nn) - 1) / nn
      pos <- data.frame(id = nodes, x = cos(ang), y = sin(ang))
      e <- cbind(kp, stats::setNames(pos[match(kp$id1, pos$id), 2:3], c("x0", "y0")), stats::setNames(pos[match(kp$id2, pos$id), 2:3], c("x1", "y1")))
      deg <- c(first = fig_label("first_degree", "First-degree", opts), second = fig_label("second_degree", "Second-degree", opts))
      e$degree <- factor(ifelse(e$r >= st$kin_thresholds[1], deg[["first"]], deg[["second"]]), levels = deg)
      # same colours as the classes they stand for: first-degree = full-sib, second-degree = half-sib
      dcol <- stats::setNames(rel_col[c(rel_lab[["fullsib"]], rel_lab[["halfsib"]])], deg)
      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = e, ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1, colour = degree, linetype = degree, linewidth = r), alpha = 0.8) +
        ggplot2::geom_point(data = pos, ggplot2::aes(x, y), size = 1.2, colour = eco_pal$ink) +
        ggplot2::scale_linewidth(range = c(0.3, 1.2), guide = "none") +
        ggplot2::scale_colour_manual(values = dcol, drop = TRUE) +
        ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed"), deg), drop = TRUE) +
        ggplot2::coord_equal(clip = "off") + ggplot2::labs(x = NULL, y = NULL, colour = NULL, linetype = NULL) +
        # positions on the circle carry no meaning: no axes (theme() is replaced by the core theme, guides are not)
        ggplot2::guides(x = "none", y = "none", colour = ggplot2::guide_legend(position = "bottom"), linetype = ggplot2::guide_legend(position = "bottom"))
      if (nn <= 30) p <- p + ggplot2::geom_text(data = pos, ggplot2::aes(1.12 * x, 1.12 * y, label = id,
                                                 hjust = ifelse(abs(x) < 0.2, 0.5, ifelse(x > 0, 0, 1))), size = 2.2, colour = eco_pal$muted)
      figs$kin_network <- with_caption(p, sprintf("Putative close kin: %d individuals (points, on a circle) joined when their estimated relatedness is at least %s. Line width is proportional to r.%s %sCheck the misclassification rates (kin_resolution) before acting on single links.",
        nn,
        if (only_first) sprintf("%.3f (first-degree relatives only; second-degree links were too many to draw)", st$kin_thresholds[1])
        else sprintf("%.3f (second-degree, dashed) or %.3f (first-degree, solid)", st$kin_thresholds[2], st$kin_thresholds[1]),
        if (top) " Only the 300 strongest links are drawn." else "",
        if (length(unique(kp$sample)) > 1) "All samples are drawn together. " else ""))
    }
  }

  # design_precision (precomputed simulation grid in design_grid.csv)
  gf <- file.path(MODULE_DIR, "design_grid.csv")
  if (file.exists(gf)) {
    g <- utils::read.csv(gf)
    mk <- if (identical(ob$marker, "msat")) "msat" else "snp"
    g <- g[g$marker == mk, ]
    if (nrow(g)) {
      unit <- if (mk == "snp") fig_label("snps", "SNPs", opts) else fig_label("microsatellites", "microsatellites", opts)
      g$loci <- factor(paste(g$n_loci, unit), levels = paste(sort(unique(g$n_loci)), unit))
      g$Ne_lab <- factor(paste(fig_label("true_Ne", "True Ne =", opts), g$Ne), levels = paste(fig_label("true_Ne", "True Ne =", opts), sort(unique(g$Ne))))
      p <- ggplot2::ggplot(g, ggplot2::aes(n_ind, cv, colour = loci, shape = loci, linetype = loci)) +
        ggplot2::geom_hline(yintercept = c(0.2, 0.5), linetype = 3, colour = eco_pal$null) +
        ggplot2::geom_line(linewidth = 0.6) + ggplot2::geom_point(size = 1.6) +
        ggplot2::facet_wrap(~Ne_lab) + ggplot2::scale_y_log10(breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2)) +
        ggplot2::scale_x_log10(breaks = sort(unique(g$n_ind)), expand = ggplot2::expansion(mult = 0.12)) +
        ggplot2::scale_colour_manual(values = pg_entity_cols(levels(g$loci), opts)) +
        ggplot2::labs(x = fig_label("n_ind", "Individuals sampled (log scale)", opts), y = fig_label("cv_inv_ne", "CV of 1/Ne (log scale)", opts),
                      colour = NULL, shape = NULL, linetype = NULL) +
        ggplot2::guides(colour = ggplot2::guide_legend(position = "bottom"), shape = ggplot2::guide_legend(position = "bottom"), linetype = ggplot2::guide_legend(position = "bottom"))
      us <- do.call(rbind, lapply(ob$ld, function(r) data.frame(n_ind = r$n_ind, n_loci = r$n_loci)))
      if (!is.null(us)) p <- p + ggplot2::geom_vline(data = unique(us["n_ind"]), ggplot2::aes(xintercept = n_ind), colour = thr_col, linewidth = 0.4, linetype = 2)
      attr(p, "fig_size") <- pg_size("double", 95, opts)
      figs$design_precision <- with_caption(p, sprintf("Precision of the LD estimate by design, from simulated Wright-Fisher populations (ecotools validation, %s): coefficient of variation of 1/Ne (the quantity LDNe estimates nearly without bias) against the number of individuals sampled, for several numbers of loci (colour, symbol, line type) and true Ne (panels). Dotted lines: CV 0.2 and 0.5. %s dashed line: this study's sample size%s.",
        if (mk == "snp") "SNPs" else "microsatellites", if (grey) "Black" else "Red",
        if (!is.null(us)) sprintf(" (%s individuals, %s loci after filtering)", paste(unique(us$n_ind), collapse = ", "), paste(unique(us$n_loci), collapse = ", ")) else ""))
    }
  }
  figs
}
