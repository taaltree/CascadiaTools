# Environment setup shared by every CascadiaTools script.
#
# The toolkit is relocatable: its root is CASCADIATOOLS_ROOT (or the older
# ECOTOOLS_ROOT), which core/run.R and core/load.R set from their own location
# before this file is sourced. Library paths are set explicitly because some
# sites' Renviron overrides R_LIBS_USER.

CASCADIA_NAME <- "CascadiaTools"
CASCADIA_URL <- "https://cascadiawildlifelab.org/tools/"
CASCADIA_CITE <- paste("Levi, T. CascadiaTools: validated analysis modules and protocols for",
                       "AI-assisted ecological analysis. Cascadia Wildlife Lab, Oregon State University.",
                       "https://github.com/taaltree/CascadiaTools")

eco_root <- function() {
  for (v in c("CASCADIATOOLS_ROOT", "ECOTOOLS_ROOT")) {
    r <- Sys.getenv(v, "")
    if (nzchar(r)) return(normalizePath(r, mustWork = FALSE))
  }
  if (file.exists(file.path(getwd(), "core", "run.R"))) return(normalizePath(getwd()))
  stop("Cannot find the CascadiaTools folder: set CASCADIATOOLS_ROOT, or run from it.", call. = FALSE)
}

eco_setup_libs <- function() {
  root <- eco_root()
  # A private library inside the toolkit (Rlib/) is used when present, for
  # example on a shared cluster; otherwise packages come from the user's
  # normal R library (installed by tools/install_packages.R).
  # CASCADIATOOLS_RLIB adds a shared library (e.g. a lab's cluster install),
  # separated by ':' (';' on Windows) for several.
  extra <- strsplit(Sys.getenv("CASCADIATOOLS_RLIB", ""), .Platform$path.sep, fixed = TRUE)[[1]]
  libs <- c(file.path(root, "Rlib"), extra)
  libs <- libs[nzchar(libs) & dir.exists(libs)]
  if (getRversion() < "4.5.0")
    warning("CascadiaTools needs R 4.5 or newer; running ", R.version.string)
  .libPaths(unique(c(libs, .libPaths())))
  invisible(.libPaths())
}

# Python for the NumPyro engines: CASCADIATOOLS_PYTHON, else a .venv inside
# the toolkit, else the parent folder's .venv (lab cluster), else python3.
eco_python <- function() {
  for (v in c("CASCADIATOOLS_PYTHON", "ECOTOOLS_PYTHON")) {
    p <- Sys.getenv(v, "")
    if (nzchar(p)) return(p)
  }
  bin <- if (.Platform$OS.type == "windows") file.path("Scripts", "python.exe") else file.path("bin", "python")
  for (d in c(file.path(eco_root(), ".venv"), file.path(dirname(eco_root()), ".venv")))
    if (file.exists(file.path(d, bin))) return(file.path(d, bin))
  "python3"
}

eco_require <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing R packages: ", paste(miss, collapse = ", "),
         ". Install them with: Rscript tools/install_packages.R", call. = FALSE)
  invisible(TRUE)
}

eco_setup_libs()
