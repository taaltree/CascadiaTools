# popgen-ne: effective population size (LD and temporal methods), pairwise
# relatedness / close kin, heterozygosity and an optional bottleneck test.
# Engine: R implementation of the NeEstimator v2 formulas (engine.R),
# cross-checked against the NeEstimator C program (VALIDATION.md).

for (f in c("engine.R", "simulate.R", "io.R", "output.R")) sys.source(file.path(MODULE_DIR, f), envir = environment())

pg_label <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# --- fit -------------------------------------------------------------------
fit <- function(data, settings) {
  st <- merge_list(list(pcrit = "auto", jackknife = TRUE, mating = "random", max_loci = 3000,
                        temporal = TRUE, relatedness = TRUE, rel_primary = "LR",
                        kin_thresholds = c(first = 0.354, second = 0.177), kin_sim_pairs = 300,
                        pcrit_grid = c(0.01, 0.02, 0.05, 0.1), heterozygosity = TRUE,
                        bottleneck = FALSE, ci = 0.95, deff_c = 0.5), settings)
  st$kin_thresholds <- sort(unlist(st$kin_thresholds), decreasing = TRUE)
  st$pcrit_grid <- unlist(st$pcrit_grid)
  ind <- data$ind; geno <- data$geno
  rows <- list(); notes <- character(); obj <- list(ld = list(), temporal = list(), rel = list(), het = list())
  add <- function(param, est, sd, lo, hi, ci_type) rows[[length(rows) + 1]] <<-
    data.frame(param = param, estimate = est, sd = sd, lower = lo, upper = hi, ci_type = ci_type, stringsAsFactors = FALSE)
  samples <- unique(ind$sample)
  one <- length(samples) == 1
  # ---- LD Ne per sample ----
  for (s in samples) {
    idx <- which(ind$sample == s & ind$use_for_ne)
    if (length(idx) < 10) { notes <- c(notes, sprintf("Sample %s: fewer than 10 individuals; no LD estimate.", s)); next }
    r <- ld_ne(geno, idx, pcrit = st$pcrit, jackknife = st$jackknife, mating = st$mating,
               max_loci = st$max_loci, ci = st$ci)
    lab <- pg_label(s)
    sd_h <- (r$inv_hi - r$inv_lo) / 2 / (2 * stats::qnorm(1 - (1 - st$ci) / 2))
    add(paste0("Ne_LD_", lab), r$ne, NA, r$ne_lo, r$ne_hi, r$ci_type)
    add(paste0("het_loss_LD_", lab), r$inv / 2, sd_h, r$inv_lo / 2, r$inv_hi / 2, r$ci_type)
    if (length(r$note)) notes <- c(notes, paste0("Sample ", s, ": ", r$note))
    r$sample <- s; r$n_ind <- length(idx)
    if (length(st$pcrit_grid)) {
      r$by_pcrit <- do.call(rbind, lapply(st$pcrit_grid, function(pc) {
        z <- tryCatch(ld_ne(geno, idx, pcrit = pc, jackknife = st$jackknife, mating = st$mating, max_loci = st$max_loci, ci = st$ci),
                      error = function(e) NULL)
        if (is.null(z) || is.na(z$inv)) return(NULL)
        data.frame(pcrit = pc, ne = z$ne, ne_lo = z$ne_lo, ne_hi = z$ne_hi, n_loci = z$n_loci)
      }))
    }
    obj$ld[[s]] <- r
    if (isTRUE(st$heterozygosity)) {
      h <- het_summary(geno, idx, st$ci)
      for (k in seq_len(nrow(h$table)))
        add(paste0(h$table$stat[k], "_", lab), h$table$estimate[k], h$table$sd[k], h$table$lower[k], h$table$upper[k], "jackknife over loci")
      if (isTRUE(st$bottleneck)) h$bottleneck <- bottleneck_test(h$per_locus)
      obj$het[[s]] <- h
    }
  }
  # ---- temporal Ne per population (earliest vs latest sample) ----
  if (isTRUE(st$temporal)) {
    pops <- unique(ind$pop)
    for (pp in pops) {
      tt_all <- table(ind$time[ind$pop == pp & ind$use_for_ne])
      tt <- sort(as.numeric(names(tt_all)[tt_all >= 10]))
      if (any(tt_all < 10)) notes <- c(notes, sprintf("Population %s: time point(s) %s have fewer than 10 individuals and were not used by the temporal method.",
                                                      pp, paste(names(tt_all)[tt_all < 10], collapse = ", ")))
      if (length(tt) < 2) next
      i1 <- which(ind$pop == pp & ind$time == tt[1] & ind$use_for_ne)
      i2 <- which(ind$pop == pp & ind$time == tt[length(tt)] & ind$use_for_ne)
      gens <- (tt[length(tt)] - tt[1]) / data$generation_time
      cn <- if (!is.null(data$census_N) && pp %in% names(data$census_N)) data$census_N[[pp]] else NA
      # drift LD for the design effect: LD estimates of the two samples (if any)
      ss <- unique(ind$sample[c(i1, i2)])
      r2d <- mean(vapply(ss, function(x) if (!is.null(obj$ld[[x]]) && is.finite(obj$ld[[x]]$rprime)) obj$ld[[x]]$rprime else NA_real_, 0), na.rm = TRUE)
      if (!is.finite(r2d)) r2d <- 0
      tm <- temporal_ne(geno, i1, i2, gens, census_N = cn, pcrit = st$pcrit, ci = st$ci, r2_drift = r2d, deff_c = st$deff_c)
      if (is.null(tm)) { notes <- c(notes, sprintf("Population %s: temporal method not possible (too few usable loci).", pp)); next }
      suf <- if (length(pops) == 1) "" else paste0("_", pg_label(pp))
      for (m in c("JR", "NT", "Pollak")) {
        h <- temporal_rows(tm, m, gens)
        nm <- if (m == "JR") "" else paste0("_", m)
        add(paste0("Ne_temporal", nm, suf), inv_to_ne(2 * h[["het"]]), NA, inv_to_ne(2 * h[["het_hi"]]), inv_to_ne(2 * h[["het_lo"]]), tm[[m]]$ci_type)
        add(paste0("het_loss_temporal", nm, suf), h[["het"]], (h[["het_hi"]] - h[["het_lo"]]) / (2 * stats::qnorm(1 - (1 - st$ci) / 2)),
            h[["het_lo"]], h[["het_hi"]], tm[[m]]$ci_type)
      }
      tm$pop <- pp; tm$times <- tt[c(1, length(tt))]
      obj$temporal[[pp]] <- tm
    }
  }
  # ---- relatedness within each sample ----
  if (isTRUE(st$relatedness)) {
    kin <- data$known_kin
    for (s in samples) {
      idx <- which(ind$sample == s)
      if (length(idx) < 5) next
      rp <- rel_prep(geno, idx)
      rm <- rel_matrix(rp)
      ids <- ind$id[idx]
      ut <- which(upper.tri(rm$LR), arr.ind = TRUE)
      pairs <- data.frame(id1 = ids[ut[, 1]], id2 = ids[ut[, 2]], LR = rm$LR[ut], QG = rm$QG[ut], stringsAsFactors = FALSE)
      pairs$r <- pairs[[st$rel_primary]]
      th <- st$kin_thresholds
      pairs$kin_class <- ifelse(pairs$r >= th[1], "first-degree (full-sib or parent-offspring)",
                                ifelse(pairs$r >= th[2], "second-degree (e.g. half-sib)", "unrelated or distant"))
      res <- list(sample = s, pairs = pairs, n_loci = rp$L)
      if (!is.null(kin)) {
        kk <- kin[kin$id1 %in% ids & kin$id2 %in% ids, , drop = FALSE]
        for (cl in unique(kk$relationship)) {
          sub <- kk[kk$relationship == cl, ]
          cm <- rel_class_mean(rp, match(sub$id1, ids), match(sub$id2, ids), st$ci)
          suf <- if (sum(vapply(samples, function(x) any(kin$id1 %in% ind$id[ind$sample == x]), TRUE)) > 1) paste0("_", pg_label(s)) else ""
          for (k in seq_len(nrow(cm)))
            add(paste0("relatedness_", pg_label(cl), "_", cm$estimator[k], suf), cm$estimate[k], cm$sd[k], cm$lower[k], cm$upper[k],
                sprintf("jackknife over loci (%d pairs)", cm$n_pairs[k]))
          res$known[[cl]] <- cm
        }
        pk <- paste(pairs$id1, pairs$id2); kk2 <- c(paste(kk$id1, kk$id2), paste(kk$id2, kk$id1))
        pairs$known <- kk$relationship[match(pk, paste(kk$id1, kk$id2))]
        pairs$known[is.na(pairs$known)] <- kk$relationship[match(pk[is.na(pairs$known)], paste(kk$id2, kk$id1))]
        res$pairs <- pairs
      }
      if (st$kin_sim_pairs > 0) {
        sp <- rel_sim_pairs(rp, st$kin_sim_pairs)
        sp$r <- sp[[st$rel_primary]]
        res$sim <- sp
        res$error <- do.call(rbind, lapply(split(sp, sp$class), function(d) data.frame(
          class = d$class[1], mean_r = mean(d$r), sd_r = stats::sd(d$r),
          p_first = mean(d$r >= th[1]), p_second = mean(d$r >= th[2] & d$r < th[1]), p_unrel = mean(d$r < th[2]))))
      }
      obj$rel[[s]] <- res
    }
  }
  est <- do.call(rbind, rows)
  if (is.null(est)) stop("no estimates could be computed; check sample sizes and loci")
  list(estimates = est,
       diagnostics = list(converged = all(is.finite(est$estimate[grepl("^het_loss", est$param)])),
                          notes = notes),
       object = c(obj, list(settings = st, marker = data$marker)))
}

# --- simulate ----------------------------------------------------------------
simulate <- function(params, design, seed) {
  sim <- wf_simulate(params, design, seed)
  truth <- c(het_loss_LD_early = 1 / (2 * params$Ne), het_loss_LD_late = 1 / (2 * params$Ne),
             het_loss_temporal = 1 / (2 * params$Ne))
  for (cl in names(sim$rel_truth)) for (e in c("LR", "QG"))
    truth[paste0("relatedness_", cl, "_", e)] <- sim$rel_truth[[cl]]
  sim$rel_truth <- NULL; sim$rel_truth_alt <- NULL
  list(data = sim, truth = truth)
}

module <- ecomodule(
  id = "popgen-ne",
  title = "Effective population size and relatedness",
  description = paste("How large is the genetically effective population (Ne), how fast is it losing genetic",
                      "diversity, and which sampled individuals are close relatives? Single-sample linkage-disequilibrium",
                      "Ne (LDNe), temporal Ne from allele-frequency change between samples taken generations apart,",
                      "pairwise relatedness (Queller-Goodnight, Lynch-Ritland) with close-kin identification, and",
                      "heterozygosity summaries from microsatellite or SNP genotypes."),
  defaults = list(
    params = list(Ne = 100),
    design = list(marker = "snp", n_loci = 1000, n_ind = 60, generations_between = 5, burn_in = 15,
                  n_alleles = 10, allele_dirichlet = 1, snp_maf_min = 0.05,
                  n_fullsib_pairs = 15, n_halfsib_pairs = 15, n_unrelated_pairs = 300, missing_rate = 0)),
  simulate = simulate, fit = fit,
  settings = list(),
  fast_settings = list(kin_sim_pairs = 0, pcrit_grid = numeric(0), heterozygosity = FALSE),
  read_data = read_genotypes, check_data = check_genotypes,
  interpret = pg_interpret, figures = pg_figures, plots = function(fit, data) pg_figures(fit, data, fig_opts()),
  methods_text = pg_methods,
  power_params = c("het_loss_LD_late", "het_loss_temporal"),
  null_values = c(het_loss_LD_late = 0, het_loss_LD_early = 0, het_loss_temporal = 0),
  engine = "R implementation of NeEstimator v2 formulas (LDNe, temporal JR/NT/Pollak); moment relatedness estimators",
  references = pg_references
)
