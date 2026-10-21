# occu-waddle: simulate two-species conditional occupancy data.
#
# params (probabilities at the average covariate value, link-scale slopes):
#   occ_dominant, occ_dominant_beta_habitat                         psiA
#   occ_sub_given_dom_present, occ_sub_given_dom_present_beta_habitat  psiBA
#   occ_sub_given_dom_absent,  occ_sub_given_dom_absent_beta_habitat   psiBa
#   det_dominant_alone, det_dominant_alone_beta_effort               pA
#   det_sub_alone, det_sub_alone_beta_effort                         pB
#   det_dominant_with_sub, det_sub_with_dom_detected, det_sub_with_dom_missed   rA, rBA, rBa
# design: n_sites, n_surveys, prop_missing (share of surveys not done),
#   det_interaction (TRUE: rA, rBA, rBa estimated; FALSE: tied to pA, pB),
#   occ_interaction (TRUE: psiBA free; FALSE: psiBA tied to psiBa in the fit).
# habitat is a standardised site covariate, effort a standardised survey covariate.

wd_simulate <- function(params, design, seed) {
  set.seed(seed)
  M <- design$n_sites; J <- design$n_surveys
  zstd <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(as.vector(x), na.rm = TRUE)
  w <- matrix(stats::rbinom(M * J, 1, 1 - design$prop_missing), M, J)
  for (i in seq_len(M)) if (sum(w[i, ]) == 0) w[i, sample.int(J, 1)] <- 1
  habitat <- zstd(stats::rnorm(M))
  effort <- matrix(stats::rnorm(M * J), M, J); effort[w == 0] <- NA; effort <- zstd(effort)
  g <- function(k, d = 0) if (is.null(params[[k]])) d else params[[k]]
  lin <- function(nm, x, cov) stats::plogis(stats::qlogis(g(nm)) + g(paste0(nm, "_beta_", cov)) * x)
  psiA <- lin("occ_dominant", habitat, "habitat")
  psiBA <- lin("occ_sub_given_dom_present", habitat, "habitat")
  psiBa <- lin("occ_sub_given_dom_absent", habitat, "habitat")
  occ_int <- !isFALSE(design$occ_interaction)
  if (!occ_int) psiBA <- psiBa
  det_int <- !isFALSE(design$det_interaction)
  e0 <- effort; e0[is.na(e0)] <- 0
  pA <- lin("det_dominant_alone", e0, "effort"); pB <- lin("det_sub_alone", e0, "effort")
  if (det_int) {
    rA <- matrix(g("det_dominant_with_sub"), M, J)
    rBA <- matrix(g("det_sub_with_dom_detected"), M, J); rBa <- matrix(g("det_sub_with_dom_missed"), M, J)
  } else { rA <- pA; rBA <- pB; rBa <- pB }
  zA <- stats::rbinom(M, 1, psiA)
  zB <- stats::rbinom(M, 1, ifelse(zA == 1, psiBA, psiBa))
  both <- matrix(zA * zB, M, J)                 # site states spread over surveys
  pa <- ifelse(both == 1, rA, pA) * matrix(zA, M, J)
  yA <- matrix(stats::rbinom(M * J, 1, pa), M, J)
  pb <- ifelse(both == 1, ifelse(yA == 1, rBA, rBa), pB) * matrix(zB, M, J)
  yB <- matrix(stats::rbinom(M * J, 1, pb), M, J)
  yA[w == 0] <- NA; yB[w == 0] <- NA
  forms <- list(psiA = "~ habitat", psiBA = if (occ_int) "~ habitat" else NULL, psiBa = "~ habitat",
                pA = "~ effort", pB = "~ effort")
  if (det_int) forms <- c(forms, list(rA = "~ 1", rBA = "~ 1", rBa = "~ 1"))
  truth <- c(occ_dominant = g("occ_dominant"), occ_dominant_beta_habitat = g("occ_dominant_beta_habitat"))
  if (occ_int) truth <- c(truth, occ_sub_given_dom_present = g("occ_sub_given_dom_present"),
                          occ_sub_given_dom_present_beta_habitat = g("occ_sub_given_dom_present_beta_habitat"))
  truth <- c(truth, occ_sub_given_dom_absent = g("occ_sub_given_dom_absent"),
             occ_sub_given_dom_absent_beta_habitat = g("occ_sub_given_dom_absent_beta_habitat"),
             det_dominant_alone = g("det_dominant_alone"), det_dominant_alone_beta_effort = g("det_dominant_alone_beta_effort"),
             det_sub_alone = g("det_sub_alone"), det_sub_alone_beta_effort = g("det_sub_alone_beta_effort"))
  if (det_int) truth <- c(truth, det_dominant_with_sub = g("det_dominant_with_sub"),
                          det_sub_with_dom_detected = g("det_sub_with_dom_detected"),
                          det_sub_with_dom_missed = g("det_sub_with_dom_missed"))
  pBA0 <- if (occ_int) g("occ_sub_given_dom_present") else g("occ_sub_given_dom_absent")
  der <- wd_derive(g("occ_dominant"), pBA0, g("occ_sub_given_dom_absent"),
                   if (det_int) g("det_sub_with_dom_detected") else 0.5,
                   if (det_int) g("det_sub_with_dom_missed") else 0.5)[1, ]
  keep <- c("occ_subordinate", "occ_both", if (occ_int) c("species_interaction_factor", "occ_log_odds_ratio_sub"),
            if (det_int) "det_log_odds_ratio_sub")
  truth <- c(truth, der[keep])
  data <- list(species = c(dominant = "dominant", subordinate = "subordinate"),
               yA = yA, yB = yB, site_covs = data.frame(habitat = habitat),
               obs_covs = list(effort = effort), formulas = forms,
               site_ids = sprintf("S%03d", seq_len(M)))
  list(data = data, truth = truth)
}
