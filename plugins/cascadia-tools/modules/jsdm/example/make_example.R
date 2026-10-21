# Builds the example dataset (SIMULATED, for demonstration only): forest
# breeding birds at 120 point-count stations in a Cascade Range landscape,
# 3 visits per station. Written in the formats field crews keep:
#   bird_detections.csv  long: one row per detection record (station, visit, date, species, count)
#   stations.csv         one row per station: coordinates (UTM m) and covariates
#   bird_traits.csv      one row per species: body mass, family, order, common name
# Truth comes from the HMSC generative model (probit, traits -> responses,
# 2 latent site factors); saved in example/truth.csv for the record.
source(file.path(Sys.getenv("ECOTOOLS_ROOT", "."), "core", "load.R"))
set.seed(20261019)
dir <- file.path(eco_root(), "modules", "jsdm", "example")
spp <- data.frame(
  common = c("Hermit Warbler", "Wilson's Warbler", "MacGillivray's Warbler", "Pacific Wren",
             "Chestnut-backed Chickadee", "Red-breasted Nuthatch", "Brown Creeper", "Golden-crowned Kinglet",
             "Hammond's Flycatcher", "Swainson's Thrush", "Varied Thrush", "Dark-eyed Junco",
             "Spotted Towhee", "Western Tanager", "Hairy Woodpecker"),
  species = c("Setophaga occidentalis", "Cardellina pusilla", "Geothlypis tolmiei", "Troglodytes pacificus",
              "Poecile rufescens", "Sitta canadensis", "Certhia americana", "Regulus satrapa",
              "Empidonax hammondii", "Catharus ustulatus", "Ixoreus naevius", "Junco hyemalis",
              "Pipilo maculatus", "Piranga ludoviciana", "Dryobates villosus"),
  family = c("Parulidae", "Parulidae", "Parulidae", "Troglodytidae", "Paridae", "Sittidae", "Certhiidae",
             "Regulidae", "Tyrannidae", "Turdidae", "Turdidae", "Passerellidae", "Passerellidae",
             "Cardinalidae", "Picidae"),
  order = c(rep("Passeriformes", 14), "Piciformes"),
  body_mass_g = c(9.3, 7.5, 10.0, 9.0, 9.8, 10.1, 8.4, 6.1, 10.2, 30.8, 77.6, 19.6, 40.6, 28.1, 66.0))
ns <- nrow(spp); n <- 120
# stations: clustered transects across a 25 x 18 km landscape
cx <- runif(12, 560000, 585000); cy <- runif(12, 4880000, 4898000)
e <- rep(cx, each = 10) + rnorm(n, 0, 700); nn <- rep(cy, each = 10) + rnorm(n, 0, 700)
elev <- round(400 + 1400 * (nn - 4880000) / 18000 + rnorm(n, 0, 150))
temp <- round(24 - 0.0055 * elev + rnorm(n, 0, 0.6), 1)            # July mean max temperature (deg C)
canopy <- round(pmin(98, pmax(5, 60 + 25 * scale(rnorm(n) + 0.3 * scale(elev))[, 1])))  # canopy cover (%)
st <- data.frame(station = sprintf("CAS-%03d", seq_len(n)), utm_e = round(e), utm_n = round(nn),
                 elevation_m = elev, july_temp_c = temp, canopy_cover_pct = canopy,
                 survey_year = 2025)
# HMSC truth (standardised covariates and trait)
zs <- function(v) (v - mean(v)) / sd(v)
X <- cbind(1, zs(temp), zs(canopy))
tz <- zs(log(spp$body_mass_g))
Gamma <- rbind(c(-0.4, 0.35), c(0.5, -0.45), c(0.4, 0.3))       # rows: intercept, temp, canopy; cols: int, log mass
B <- Gamma %*% rbind(1, tz) + matrix(rnorm(3 * ns, 0, 0.5), 3, ns)
Lam <- rbind(rnorm(ns, 0, 0.9), rnorm(ns, 0, 0.6))
Eta <- matrix(rnorm(n * 2), n, 2)
L <- X %*% B + Eta %*% Lam
Z <- (L + matrix(rnorm(n * ns), n, ns) > 0) * 1
# imperfect detection over 3 visits; counts of individuals on detection
rows <- list()
for (i in seq_len(n)) for (v in 1:3) {
  date <- format(as.Date("2025-05-20") + (v - 1) * 14 + sample(0:3, 1))
  for (j in which(Z[i, ] == 1)) if (runif(1) < 0.7)
    rows[[length(rows) + 1]] <- data.frame(station = st$station[i], visit = v, date = date,
                                           species = spp$species[j], count = 1 + rpois(1, 0.6))
}
det <- do.call(rbind, rows)
utils::write.csv(st, file.path(dir, "stations.csv"), row.names = FALSE)
utils::write.csv(det, file.path(dir, "bird_detections.csv"), row.names = FALSE)
utils::write.csv(spp[, c("species", "common", "family", "order", "body_mass_g")], file.path(dir, "bird_traits.csv"), row.names = FALSE)
tv <- data.frame(param = c(paste0("beta_", rep(c("intercept", "july_temp_c", "canopy_cover_pct"), each = ns), "_", rep(spp$species, 3)),
                           "comm_mean_july_temp_c", "comm_mean_canopy_cover_pct", "trait_log_body_mass_g_on_july_temp_c"),
                 truth = c(as.vector(t(B)), Gamma[2, 1], Gamma[3, 1], Gamma[2, 2]))
utils::write.csv(tv, file.path(dir, "truth.csv"), row.names = FALSE)
cat("stations", n, "detections", nrow(det), "naive prevalence:\n"); print(round(colMeans(Z), 2))
