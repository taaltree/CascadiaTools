# Readers for the genotype formats researchers actually have.
#
# config$data keys (paths relative to the config file):
#   genotypes:  CSV/TSV genotype table, one row per individual (or use vcf:)
#   vcf:        VCF file (vcfR); populations then come from `metadata`
#   format:     auto | two_column | slash | snp012     (default auto)
#   id_col:     individual id column (default "individual")
#   pop_col:    population column (default "population"); can also be in metadata
#   loci:       optional list of locus columns (two_column: the FIRST column of each pair)
#   exclude:    optional columns to ignore (e.g. sex, date)
#   missing:    codes meaning missing (default 0, -9, -1, NA, "", "?", "./.", "N")
#   marker:     microsat | snp (default: snp for snp012/VCF, else microsat)
#   metadata:   CSV with coordinates: per population (pop_col + coords) or per
#               individual (id_col [+ pop_col] + coords; averaged per population)
#   lat_col, lon_col: decimal-degree columns (default latitude/longitude, lat/lon)
#   x_col, y_col: projected coordinates IN KILOMETRES (alternative to lat/lon)
#   ibd_distance: log (2-D habitat, default) | linear (1-D habitat: rivers, coastlines)

ps_missing_default <- c("0", "-9", "-1", "NA", "", "?", "./.", ".|.", ".", "N", "000", "0/0")

ps_read_table <- function(f) {
  sep <- if (grepl("\\.(tsv|txt)$", f, ignore.case = TRUE)) "\t" else ","
  utils::read.table(f, header = TRUE, sep = sep, colClasses = "character", check.names = FALSE,
                    na.strings = character(), quote = "\"", comment.char = "", strip.white = TRUE)
}

ps_parse_slash <- function(v, miss) {
  v <- trimws(v)
  na <- v %in% miss | is.na(v)
  parts <- strsplit(ifelse(na, "", v), "[/|:]")
  len <- lengths(parts)
  bad <- !na & len != 2
  a <- vapply(parts, function(z) if (length(z) == 2) trimws(z[1]) else NA_character_, "")
  b <- vapply(parts, function(z) if (length(z) == 2) trimws(z[2]) else NA_character_, "")
  half <- !na & !bad & (a %in% miss | b %in% miss)   # e.g. "152/0": one allele unscored
  a[na | bad | half] <- NA; b[na | bad | half] <- NA
  ai <- suppressWarnings(as.integer(a)); bi <- suppressWarnings(as.integer(b))
  if (any(!is.na(a) & is.na(ai)) || any(!is.na(b) & is.na(bi))) {
    # non-numeric allele names (e.g. nucleotides A/T): code them 1..k
    lev <- sort(unique(stats::na.omit(c(a, b))))
    ai <- match(a, lev); bi <- match(b, lev)
  }
  list(a = ai, b = bi, bad = which(bad), half = sum(half))
}

ps_resolve <- function(cfg, key) {
  f <- cfg$data[[key]]
  if (is.null(f)) return(NULL)
  if (!grepl("^/", f) && !is.null(cfg$.dir)) f <- file.path(cfg$.dir, f)
  if (!file.exists(f)) stop("File not found for data$", key, ": ", f)
  f
}

ps_read_data <- function(config) {
  dc <- config$data
  idc <- dc$id_col %||% "individual"; popc <- dc$pop_col %||% "population"
  miss <- c(ps_missing_default, unlist(dc$missing))
  problems <- character(); warns <- character()
  meta <- if (!is.null(dc$metadata)) ps_read_table(ps_resolve(config, "metadata")) else NULL
  vcf <- ps_resolve(config, "vcf")
  if (!is.null(vcf)) {
    eco_require("vcfR")
    v <- vcfR::read.vcfR(vcf, verbose = FALSE)
    gt <- vcfR::extract.gt(v, element = "GT")
    ids <- colnames(gt)
    loci <- make.unique(ifelse(is.na(vcfR::getID(v)) | vcfR::getID(v) == ".",
                               paste0(vcfR::getCHROM(v), "_", vcfR::getPOS(v)), vcfR::getID(v)))
    A1 <- A2 <- matrix(NA_integer_, length(ids), length(loci))
    for (l in seq_along(loci)) {
      pz <- ps_parse_slash(gt[l, ], miss = c(".", "./.", ".|.", NA))
      A1[, l] <- pz$a; A2[, l] <- pz$b
    }
    if (any(A1 > 1 | A2 > 1, na.rm = TRUE)) warns <- c(warns, "Some VCF sites are multiallelic; they are analysed as multiallelic loci.")
    colnames(A1) <- colnames(A2) <- loci
    if (is.null(meta) || !idc %in% names(meta) || !popc %in% names(meta))
      stop("With a VCF, data$metadata must be a CSV with columns '", idc, "' and '", popc, "' giving each sample's population.")
    pop <- meta[[popc]][match(ids, meta[[idc]])]
    if (anyNA(pop)) problems <- c(problems, sprintf("%d VCF samples are not in the metadata file (e.g. %s). Add them or remove them from the VCF.",
                                                    sum(is.na(pop)), paste(utils::head(ids[is.na(pop)], 3), collapse = ", ")))
    marker <- dc$marker %||% "snp"
  } else {
    gf <- ps_resolve(config, "genotypes")
    if (is.null(gf)) stop("config data$genotypes (a CSV) or data$vcf is required")
    tb <- ps_read_table(gf)
    if (!idc %in% names(tb)) stop("Genotype file has no individual column '", idc, "'. Columns are: ", paste(utils::head(names(tb), 8), collapse = ", "), " ... Set data$id_col.")
    ids <- tb[[idc]]
    if (popc %in% names(tb)) pop <- tb[[popc]]
    else if (!is.null(meta) && all(c(idc, popc) %in% names(meta))) pop <- meta[[popc]][match(ids, meta[[idc]])]
    else stop("No population column '", popc, "' in the genotype file or metadata. Genotype file columns are: ",
              paste(utils::head(names(tb), 8), collapse = ", "), if (ncol(tb) > 8) sprintf(" (and %d more)", ncol(tb) - 8) else "",
              ". Set data$pop_col to the column that holds the population / site name.")
    cols <- setdiff(names(tb), c(idc, popc, unlist(dc$exclude)))
    cols <- cols[nzchar(cols)]
    fmt <- dc$format %||% "auto"
    vals <- unlist(tb[cols]); vals <- vals[!vals %in% miss]
    if (fmt == "auto") {
      fmt <- if (length(vals) && mean(grepl("[/|]", vals)) > 0.9) "slash"
             else if (length(vals) && all(vals %in% c("0", "1", "2"))) "snp012" else "two_column"
      # 0 is also a common missing code: snp012 needs the user to say so
    }
    if (fmt == "slash") {
      loci <- if (!is.null(dc$loci)) unlist(dc$loci) else cols
      A1 <- A2 <- matrix(NA_integer_, length(ids), length(loci))
      for (l in seq_along(loci)) {
        pz <- ps_parse_slash(tb[[loci[l]]], miss)
        A1[, l] <- pz$a; A2[, l] <- pz$b
        if (length(pz$bad)) problems <- c(problems, sprintf("Locus %s: %d genotype(s) not of the form 'allele/allele' (e.g. individual %s has '%s').",
                                                            loci[l], length(pz$bad), ids[pz$bad[1]], tb[[loci[l]]][pz$bad[1]]))
      }
      marker <- dc$marker %||% "microsat"
    } else if (fmt == "snp012") {
      loci <- if (!is.null(dc$loci)) unlist(dc$loci) else cols
      miss012 <- setdiff(miss, "0")
      A1 <- A2 <- matrix(NA_integer_, length(ids), length(loci))
      for (l in seq_along(loci)) {
        v <- trimws(tb[[loci[l]]]); na <- v %in% miss012
        z <- suppressWarnings(as.integer(v)); z[na] <- NA
        bad <- which(!na & !(z %in% 0:2))
        if (length(bad)) problems <- c(problems, sprintf("SNP %s: value '%s' (individual %s) is not 0, 1 or 2.", loci[l], v[bad[1]], ids[bad[1]]))
        z[!(z %in% 0:2)] <- NA
        A1[, l] <- ifelse(z >= 1, 1L, 0L); A2[, l] <- ifelse(z == 2, 1L, 0L)
      }
      marker <- dc$marker %||% "snp"
    } else {
      # two columns per locus: either explicit dc$loci (first column of each pair) or consecutive pairs
      if (length(cols) %% 2) stop("Two-column format needs an even number of allele columns; found ", length(cols),
                                  ". List non-genotype columns under data$exclude.")
      first <- if (!is.null(dc$loci)) match(unlist(dc$loci), cols) else seq(1, length(cols), by = 2)
      loci <- sub("[._-](a|b|1|2|A|B)$", "", cols[first])
      loci <- make.unique(ifelse(nzchar(loci), loci, cols[first]))
      A1 <- A2 <- matrix(NA_integer_, length(ids), length(loci))
      for (l in seq_along(loci)) {
        v1 <- trimws(tb[[cols[first[l]]]]); v2 <- trimws(tb[[cols[first[l] + 1]]])
        z1 <- suppressWarnings(as.integer(v1)); z2 <- suppressWarnings(as.integer(v2))
        na <- v1 %in% miss | v2 %in% miss
        bad <- which(!na & (is.na(z1) | is.na(z2)))
        if (length(bad)) problems <- c(problems, sprintf("Locus %s: allele '%s/%s' (individual %s) is not a whole number.", loci[l], v1[bad[1]], v2[bad[1]], ids[bad[1]]))
        half <- which(xor(v1 %in% miss, v2 %in% miss))
        if (length(half)) warns <- c(warns, sprintf("Locus %s: %d individual(s) have only one allele scored; treated as missing.", loci[l], length(half)))
        z1[na | is.na(z2)] <- NA; z2[na | is.na(z1)] <- NA
        A1[, l] <- z1; A2[, l] <- z2
      }
      marker <- dc$marker %||% "microsat"
    }
    colnames(A1) <- colnames(A2) <- loci
  }
  # coordinates per population
  coords <- NULL; ctype <- "km"
  if (!is.null(meta)) {
    nm <- tolower(names(meta))
    pick <- function(opts, key) { k <- dc[[key]]; if (!is.null(k)) return(k); w <- names(meta)[nm %in% opts]; if (length(w)) w[1] else NULL }
    latc <- pick(c("latitude", "lat", "decimallatitude"), "lat_col"); lonc <- pick(c("longitude", "lon", "long", "decimallongitude"), "lon_col")
    xc <- pick(c("x_km", "x", "easting_km"), "x_col"); yc <- pick(c("y_km", "y", "northing_km"), "y_col")
    if (!is.null(latc) && !is.null(lonc)) { cx <- lonc; cy <- latc; ctype <- "lonlat" } else if (!is.null(xc) && !is.null(yc)) { cx <- xc; cy <- yc } else cx <- NULL
    if (!is.null(cx)) {
      key <- if (popc %in% names(meta) && !(idc %in% names(meta))) meta[[popc]]
             else if (idc %in% names(meta)) pop[match(meta[[idc]], ids)] else NULL
      if (is.null(key)) stop("metadata needs a '", popc, "' or '", idc, "' column to attach coordinates.")
      X <- suppressWarnings(as.numeric(meta[[cx]])); Y <- suppressWarnings(as.numeric(meta[[cy]]))
      agg <- stats::aggregate(cbind(x = X, y = Y), list(pop = key), mean, na.rm = TRUE)
      coords <- agg
      if (ctype == "lonlat" && any(abs(coords$y) > 90 | abs(coords$x) > 180, na.rm = TRUE))
        problems <- c(problems, "Latitude/longitude values out of range: are the columns swapped or projected (metres)? Use x_col/y_col in km for projected coordinates.")
      nocoord <- setdiff(unique(pop), coords$pop[is.finite(coords$x) & is.finite(coords$y)])
      unused <- setdiff(coords$pop, unique(pop))
      if (length(nocoord)) problems <- c(problems, sprintf("No coordinates for population(s): %s. Add them to the metadata file or remove data$metadata to skip isolation by distance.%s",
                                                           paste(nocoord, collapse = ", "),
                                                           if (length(unused)) sprintf(" The metadata file has population(s) that match no genotyped population: %s (a misspelling or different capitalisation?).",
                                                                                       paste(unused, collapse = ", ")) else ""))
    }
  }
  d <- ps_make_data(A1, A2, ids, pop, marker = marker, coords = coords, coord_type = ctype,
                    ibd_distance = dc$ibd_distance %||% "log",
                    source = basename(vcf %||% ps_resolve(config, "genotypes")), warnings = warns, problems = problems)
  d$warnings <- c(d$warnings, ps_data_warnings(d))
  d
}

# Non-fatal data issues (reported in the diagnostics notes).
ps_data_warnings <- function(d) {
  w <- character()
  mloc <- colMeans(is.na(d$a1))
  if (any(mloc > 0.2)) w <- c(w, sprintf("Loci with > 20%% missing genotypes: %s.", paste(d$loci[mloc > 0.2], collapse = ", ")))
  mind <- rowMeans(is.na(d$a1))
  if (any(mind > 0.5)) w <- c(w, sprintf("%d individual(s) missing > 50%% of loci (e.g. %s); consider removing them.", sum(mind > 0.5), paste(utils::head(d$ind[mind > 0.5], 3), collapse = ", ")))
  np <- table(d$pop)
  if (any(np < 10)) w <- c(w, sprintf("Populations with fewer than 10 individuals (%s): their He, FIS and pairwise FST are imprecise.", paste(names(np)[np < 10], collapse = ", ")))
  mono <- vapply(seq_along(d$loci), function(l) length(unique(stats::na.omit(c(d$a1[, l], d$a2[, l])))) < 2, TRUE)
  if (any(mono)) w <- c(w, sprintf("%d monomorphic locus/loci (%s).", sum(mono), paste(utils::head(d$loci[mono], 5), collapse = ", ")))
  w
}

ps_check_data <- function(d) {
  pr <- d$problems
  if (nlevels(d$pop) < 2) pr <- c(pr, "Only one population: FST and isolation by distance need at least two (three for IBD).")
  np <- table(d$pop); small <- names(np)[np < 3]
  if (length(small)) pr <- c(pr, sprintf("Population(s) with fewer than 3 individuals: %s. Merge them with a neighbouring sample or drop them.", paste(small, collapse = ", ")))
  pl <- levels(d$pop); key <- tolower(gsub("[[:space:]_.-]", "", pl))
  near <- split(pl, key); near <- near[lengths(near) > 1]
  if (length(near)) pr <- c(pr, sprintf("Population names that differ only in capitalisation, spaces or punctuation: %s. They are treated as different populations; correct the spelling in the genotype file.",
                                        paste(vapply(near, function(z) paste(sprintf("'%s'", z), collapse = " vs "), ""), collapse = "; ")))
  dup <- unique(d$ind[duplicated(d$ind)])
  if (length(dup)) pr <- c(pr, sprintf("Duplicated individual ids: %s. Each row must be one individual.", paste(utils::head(dup, 5), collapse = ", ")))
  allmiss <- d$loci[colSums(!is.na(d$a1)) == 0]
  if (length(allmiss)) pr <- c(pr, sprintf("Loci with no genotypes at all: %s. Remove them.", paste(allmiss, collapse = ", ")))
  typed <- vapply(seq_along(d$loci), function(l) tapply(!is.na(d$a1[, l]), d$pop, sum), numeric(nlevels(d$pop)))
  typed <- matrix(typed, nlevels(d$pop))
  if (any(typed < 2)) {
    w <- which(typed < 2, arr.ind = TRUE)
    pr <- c(pr, sprintf("Population %s has fewer than 2 typed individuals at locus %s (and %d other population x locus combinations). Remove that locus or population.",
                        levels(d$pop)[w[1, 1]], d$loci[w[1, 2]], nrow(w) - 1))
  }
  pr
}
