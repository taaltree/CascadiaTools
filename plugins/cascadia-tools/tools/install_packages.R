# Install the R packages CascadiaTools needs into your normal R library.
#
#   Rscript tools/install_packages.R            # everything (about 20-40 min the first time)
#   Rscript tools/install_packages.R --check    # only report what is missing
#
# Nimble-based modules compile C++, so you also need a compiler:
#   Windows: Rtools (https://cran.r-project.org/bin/windows/Rtools/)
#   macOS:   Xcode command-line tools (run `xcode-select --install` in Terminal)
#   Linux:   build-essential / gcc-c++ from your package manager

pkgs <- c(
  # core and figures
  "yaml", "jsonlite", "ggplot2", "patchwork", "scales", "ggrepel", "ragg", "svglite", "base64enc",
  "posterior", "coda", "rlang",
  # engines
  "nimble", "nimbleEcology", "unmarked", "spOccupancy", "secr", "openCR", "glmmTMB", "lme4", "DHARMa",
  "performance", "marginaleffects", "MASS", "sandwich", "reformulas", "AICcmodavg",
  "dagitty", "ggraph", "igraph", "lavaan", "piecewiseSEM",
  "momentuHMM", "crawl", "amt", "ctmm", "sf", "terra", "geosphere", "lubridate", "survival", "glmnet",
  "marked", "R2ucare", "Distance", "mrds", "IPMbook", "MARSS", "popbio", "maxnet",
  "vegan", "permute", "betapart", "indicspecies", "Hmsc", "ape", "ade4",
  "adegenet", "hierfstat", "pegas", "poppr", "vcfR", "FSA"
)

check_only <- "--check" %in% commandArgs(trailingOnly = TRUE)
if (getRversion() < "4.5.0")
  stop("CascadiaTools needs R 4.5 or newer (you have ", R.version.string,
       "). Download R from https://cran.r-project.org", call. = FALSE)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
cat(sprintf("%d of %d packages already installed.\n", length(pkgs) - length(missing), length(pkgs)))
if (!length(missing)) { cat("Everything is installed.\n"); quit(status = 0) }
cat("Missing:", paste(missing, collapse = ", "), "\n")
if (check_only) quit(status = 1)

has_compiler <- nzchar(Sys.which(if (.Platform$OS.type == "windows") "make" else "g++")) ||
  nzchar(Sys.which("clang++"))
if (!has_compiler)
  message("No C++ compiler found. Packages will install, but Nimble models will not compile ",
          "until you install one (see the top of this file).")

options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(missing, Ncpus = max(1L, parallel::detectCores() - 1L))
still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
if (length(still)) {
  cat("\nCould not install:", paste(still, collapse = ", "), "\n",
      "Modules that need them will say so; everything else works.\n")
  quit(status = 1)
}
cat("\nAll packages installed. Try: Rscript core/run.R list\n")
