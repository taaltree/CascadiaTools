# Cross-check: Hmsc (module fit) vs an independent likelihood JSDM in glmmTMB
# (reduced-rank latent factors, rr(); McGillycuddy et al. 2025) on the same
# simulated datasets. Probit presence/absence, 2 latent factors.
#   scripts/ejob -c 16 -m 32 -- scripts/Rscript modules/jsdm/validation/crosscheck_glmmtmb.R [nrep] [label]
source(file.path(Sys.getenv("ECOTOOLS_ROOT", "."), "core", "load.R"))
m <- load_module("jsdm")
eco_require(c("glmmTMB"))
pair_names <- get("pair_names", envir = environment(m$simulate))
a <- commandArgs(TRUE)
nrep <- if (length(a)) as.integer(a[1]) else 16
lab <- if (length(a) > 1) a[2] else "default"
outdir <- file.path(eco_root(), "results", "jsdm", "crosscheck"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tmb_jsdm <- function(dat, d = 2) {
  Y <- dat$Y; sp <- colnames(Y); n <- nrow(Y); ns <- ncol(Y)
  sd <- dat$site_data
  z <- function(v) (v - mean(v)) / sd(v)
  long <- data.frame(site = factor(rep(rownames(Y), ns), levels = rownames(Y)),
                     species = factor(rep(sp, each = n), levels = sp),
                     y = as.vector(Y),
                     temperature = rep(z(sd$temperature), ns), forest = rep(z(sd$forest), ns))
  f <- glmmTMB::glmmTMB(y ~ 0 + species + species:temperature + species:forest + rr(0 + species | site, d = d),
                        family = stats::binomial(link = "probit"), data = long,
                        control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000)))
  cf <- glmmTMB::fixef(f)$cond; se <- sqrt(diag(stats::vcov(f)$cond))
  nm <- names(cf)
  par <- ifelse(grepl(":", nm), paste0("beta_", sub(".*:", "", nm), "_", sub("^species([^:]*):.*", "\\1", nm)),
                paste0("beta_intercept_", sub("^species", "", nm)))
  S <- glmmTMB::VarCorr(f)$cond$site[seq_len(ns), seq_len(ns)]
  lc <- stats::cov2cor(S + diag(1, ns))
  pn <- pair_names(sp)
  # population-averaged probit slopes: beta / sqrt(1 + Omega_jj)
  slope <- grepl(":", nm)
  spn <- sub("^species([^:]*):.*", "\\1", nm[slope])
  pa <- cf[slope] / sqrt(1 + diag(S)[match(spn, sp)])
  data.frame(param = c(par, sub("^beta_", "popavg_", par[slope]), paste0("latcor_", pn)),
             tmb = c(cf, pa, lc[upper.tri(lc)]),
             tmb_se = c(se, rep(NA, sum(slope) + length(pn))),
             tmb_conv = f$fit$convergence == 0 && isTRUE(f$sdr$pdHess))
}

res <- parallel::mclapply(seq_len(nrep), function(r) {
  s <- m$simulate(m$defaults$params, m$defaults$design, 5000 + r)
  t0 <- proc.time()[[3]]
  h <- m$fit(s$data, m$fast_settings)
  th <- proc.time()[[3]] - t0
  t0 <- proc.time()[[3]]
  g <- tryCatch(tmb_jsdm(s$data), error = function(e) { message("glmmTMB failed: ", conditionMessage(e)); NULL })
  tg <- proc.time()[[3]] - t0
  if (is.null(g)) return(NULL)
  e <- merge(h$estimates[, c("param", "estimate", "sd", "lower", "upper")], g, by = "param")
  e$truth <- s$truth[e$param]; e$rep <- r; e$t_hmsc <- th; e$t_tmb <- tg
  e$hmsc_rhat <- h$diagnostics$rhat_max
  e
}, mc.cores = min(16, nrep), mc.preschedule = FALSE)
r <- do.call(rbind, Filter(Negate(is.null), res))
if (is.null(r) || !nrow(r)) stop("no successful cross-check replicates")
utils::write.csv(r, file.path(outdir, paste0("crosscheck_", lab, ".csv")), row.names = FALSE)
r$cls <- sub("_.*", "", r$param)
r$tmb_cover <- with(r, ifelse(is.na(tmb_se), NA, tmb - 1.96 * tmb_se <= truth & tmb + 1.96 * tmb_se >= truth))
r$hmsc_cover <- with(r, lower <= truth & upper >= truth)
cat(sprintf("reps: %d; glmmTMB converged (pdHess) in %d; mean time Hmsc %.0f s, glmmTMB %.0f s\n",
            length(unique(r$rep)), sum(tapply(r$tmb_conv, r$rep, all)), mean(r$t_hmsc), mean(r$t_tmb)))
for (k in c("beta", "popavg", "latcor")) {
  d <- r[r$cls == k, ]
  cat(sprintf("%-6s n=%d  cor(Hmsc, glmmTMB)=%.3f  mean|diff|=%.3f  median|diff|/postSD=%.2f  RMSE vs truth: Hmsc %.3f, glmmTMB %.3f  sign agreement %.2f  coverage Hmsc %.2f glmmTMB(Wald) %s\n",
              k, nrow(d), stats::cor(d$estimate, d$tmb), mean(abs(d$estimate - d$tmb)),
              stats::median(abs(d$estimate - d$tmb) / d$sd),
              sqrt(mean((d$estimate - d$truth)^2)), sqrt(mean((d$tmb - d$truth)^2)),
              mean(sign(d$estimate) == sign(d$tmb)), mean(d$hmsc_cover),
              if (all(is.na(d$tmb_cover))) "n/a" else sprintf("%.2f", mean(d$tmb_cover, na.rm = TRUE))))
}
