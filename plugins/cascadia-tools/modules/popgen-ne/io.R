# Reading researchers' genotype files into the module's data structure:
#   list(geno = int array [n, L, 2] (NA missing), ind = data.frame(id, pop,
#        sample, time, use_for_ne), known_kin = data.frame(id1, id2,
#        relationship) or NULL, census_N = named numeric (per pop) or NULL,
#        generation_time = number, marker = "msat"/"snp")
#
# Supported formats (data$format):
#   two_column : one row per individual, two columns per locus (GenAlEx /
#                CERVUS / COANCESTRY style), e.g. Sfo18_a, Sfo18_b
#   one_column : one column per locus holding "152/160" or "A/G" (allele_sep)
#   dosage     : SNP 0/1/2 allele counts (one column per SNP)
#   genepop    : GENEPOP .gen file (2- or 3-digit alleles), pop from popmap
#                or from the order of 'Pop' blocks
#   vcf        : VCF via vcfR (biallelic SNPs), pop from popmap

pg_resolve <- function(config, f) {
  if (is.null(f)) return(NULL)
  if (grepl("^/", f)) f else file.path(config$.dir, f)
}

pg_missing <- function(x, codes) {
  x <- trimws(as.character(x))
  x[x %in% as.character(codes) | x == "" | toupper(x) == "NA"] <- NA
  x
}

# Map character allele labels to integer codes per locus (keeps numeric
# fragment sizes as they are).
pg_codes <- function(a1, a2) {
  v <- c(a1, a2)
  if (all(is.na(v) | grepl("^-?[0-9]+$", v))) return(list(as.integer(a1), as.integer(a2)))
  lev <- sort(unique(v[!is.na(v)]))
  list(match(a1, lev), match(a2, lev))
}

pg_meta <- function(df, dd) {
  idc <- if (is.null(dd$id_column)) names(df)[1] else dd$id_column
  if (!idc %in% names(df)) stop("id_column '", idc, "' not found in the genotype file. Columns are: ",
                                paste(utils::head(names(df), 8), collapse = ", "), if (ncol(df) > 8) ", ..." else "",
                                ". Column names are case-sensitive.")
  popc <- dd$pop_column
  if (!is.null(popc) && !popc %in% names(df)) stop("pop_column '", popc, "' not found in the genotype file")
  tc <- dd$time_column; sc <- dd$sample_column; uc <- dd$use_for_ne_column
  for (cc in c(tc, sc, uc)) if (!cc %in% names(df)) stop("column '", cc, "' not found in the genotype file (column names are case-sensitive; columns: ",
                                                        paste(utils::head(names(df), 8), collapse = ", "), if (ncol(df) > 8) ", ..." else "", ")")
  # 2. sampling times must be numbers (was: silently NA with a coercion warning)
  if (!is.null(tc)) {
    tv <- trimws(as.character(df[[tc]]))
    bad <- !is.na(tv) & nzchar(tv) & is.na(suppressWarnings(as.numeric(tv)))
    if (any(bad)) stop("time_column '", tc, "' has values that are not numbers: ",
                       paste(utils::head(sprintf("'%s' (individual %s)", tv[bad], df[[idc]][bad]), 5), collapse = ", "),
                       ". Use the sampling year or generation as a number.")
  }
  ind <- data.frame(id = as.character(df[[idc]]),
                    pop = if (is.null(popc)) "pop1" else as.character(df[[popc]]),
                    time = if (is.null(tc)) NA_real_ else as.numeric(df[[tc]]), stringsAsFactors = FALSE)
  ind$use_for_ne <- if (is.null(uc)) TRUE else as.logical(toupper(as.character(df[[uc]])) %in% c("TRUE", "T", "1", "YES", "Y"))
  ind$sample <- if (!is.null(sc)) as.character(df[[sc]]) else
    if (all(is.na(ind$time))) ind$pop else ifelse(is.na(ind$time), ind$pop, paste0(ind$pop, "_", ind$time))
  list(ind = ind, meta_cols = c(idc, popc, tc, sc, uc))
}

# Stop if a locus is coded with numbers (allele sizes) but a few cells are not
# numbers (typos, notes, stray text): those would otherwise become extra alleles.
pg_check_numeric <- function(df, c1, c2, miss, ids) {
  v <- c(pg_missing(df[[c1]], miss), pg_missing(df[[c2]], miss))
  num <- grepl("^-?[0-9]+$", v)
  bad <- !is.na(v) & !num
  if (any(bad) && mean(num[!is.na(v)]) > 0.5) {
    w <- which(bad); i <- ((w - 1) %% nrow(df)) + 1
    stop("locus columns ", c1, "/", c2, ": ", sum(bad), " value(s) are not allele sizes, e.g. '", v[w[1]],
         "' for individual ", ids[i[1]], ". Correct them or code them as missing (", paste(miss[1], collapse = ""), ").")
  }
}

pg_read_table <- function(f) utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")

read_genotypes <- function(config) {
  dd <- config$data
  if (is.null(dd$genotypes)) stop("config data: needs 'genotypes' (path to the genotype file)")
  f <- pg_resolve(config, dd$genotypes)
  if (!file.exists(f)) stop("genotype file not found: ", f)
  fmt <- if (is.null(dd$format)) "two_column" else dd$format
  miss <- if (is.null(dd$missing)) c("0", "-9", "NA", "?", "-") else as.character(unlist(dd$missing))
  if (fmt %in% c("two_column", "one_column", "dosage")) {
    df <- pg_read_table(f)
    mt <- pg_meta(df, dd)
    lc <- setdiff(names(df), mt$meta_cols)
    if (!is.null(dd$first_locus_column)) {
      k <- if (is.numeric(dd$first_locus_column)) dd$first_locus_column else match(dd$first_locus_column, names(df))
      if (is.na(k) || k < 1 || k > ncol(df)) stop("first_locus_column '", dd$first_locus_column, "' not found in the genotype file (column names are case-sensitive). The columns after the metadata are: ",
                                                  paste(utils::head(setdiff(names(df), mt$meta_cols), 4), collapse = ", "), ", ...")
      lc <- names(df)[k:ncol(df)]
    }
    n <- nrow(df)
    if (fmt == "two_column") {
      if (length(lc) %% 2) stop("two_column format: found ", length(lc), " locus columns, which is odd. Each locus needs two allele columns (e.g. Loc1_a, Loc1_b). Check for an extra metadata column and list it as pop_column/time_column, or set first_locus_column.")
      L <- length(lc) / 2
      loci <- sub("([._ -]?(a|b|A|B|1|2))$", "", lc[seq(1, 2 * L, 2)])
      geno <- array(NA_integer_, c(n, L, 2))
      for (j in seq_len(L)) {
        pg_check_numeric(df, lc[2 * j - 1], lc[2 * j], miss, df[[mt$meta_cols[1]]])
        cc <- pg_codes(pg_missing(df[[lc[2 * j - 1]]], miss), pg_missing(df[[lc[2 * j]]], miss))
        geno[, j, 1] <- cc[[1]]; geno[, j, 2] <- cc[[2]]
      }
    } else if (fmt == "one_column") {
      sep <- if (is.null(dd$allele_sep)) "/" else dd$allele_sep
      L <- length(lc); loci <- lc
      geno <- array(NA_integer_, c(n, L, 2))
      for (j in seq_len(L)) {
        v <- pg_missing(df[[lc[j]]], miss)
        sp <- if (nzchar(sep)) strsplit(v, sep, fixed = TRUE) else lapply(v, function(s) if (is.na(s)) NA else strsplit(s, "")[[1]])
        a1 <- vapply(sp, function(s) if (length(s) == 2) s[1] else NA_character_, "")
        a2 <- vapply(sp, function(s) if (length(s) == 2) s[2] else NA_character_, "")
        a1 <- pg_missing(a1, miss); a2 <- pg_missing(a2, miss)
        a1[is.na(a2)] <- NA; a2[is.na(a1)] <- NA
        cc <- pg_codes(a1, a2); geno[, j, 1] <- cc[[1]]; geno[, j, 2] <- cc[[2]]
      }
    } else {
      L <- length(lc); loci <- lc
      geno <- array(NA_integer_, c(n, L, 2))
      for (j in seq_len(L)) {
        v <- suppressWarnings(as.integer(pg_missing(df[[lc[j]]], setdiff(miss, "0"))))
        if (any(!is.na(v) & !v %in% 0:2)) stop("dosage format: locus ", lc[j], " has values other than 0, 1, 2 (", paste(utils::head(unique(v[!v %in% 0:2]), 3), collapse = ", "), ")")
        geno[, j, 1] <- ifelse(is.na(v), NA, ifelse(v >= 1, 2L, 1L)); geno[, j, 2] <- ifelse(is.na(v), NA, ifelse(v == 2, 2L, 1L))
      }
    }
    ind <- mt$ind
  } else if (fmt == "genepop") {
    x <- readLines(f, warn = FALSE)
    pl <- grep("^\\s*pop\\s*$", x, ignore.case = TRUE)
    if (!length(pl)) stop("genepop: no 'Pop' lines found")
    hdr <- x[2:(pl[1] - 1)]
    loci <- trimws(unlist(strsplit(hdr, ",")))
    loci <- loci[nzchar(loci)]
    rows <- list(); popn <- integer()
    for (b in seq_along(pl)) {
      end <- if (b < length(pl)) pl[b + 1] - 1 else length(x)
      if (end <= pl[b]) next
      for (ln in x[(pl[b] + 1):end]) if (nzchar(trimws(ln))) { rows[[length(rows) + 1]] <- ln; popn <- c(popn, b) }
    }
    ids <- trimws(sub(",.*$", "", unlist(rows)))
    gts <- lapply(rows, function(r) strsplit(trimws(sub("^[^,]*,", "", r)), "\\s+")[[1]])
    L <- length(loci); n <- length(rows)
    if (any(lengths(gts) != L)) stop("genepop: some individuals do not have ", L, " genotypes")
    geno <- array(NA_integer_, c(n, L, 2))
    for (i in seq_len(n)) {
      g <- gts[[i]]; h <- nchar(g) / 2
      a1 <- as.integer(substr(g, 1, h)); a2 <- as.integer(substr(g, h + 1, 2 * h))
      a1[a1 == 0 | a2 == 0] <- NA; a2[is.na(a1)] <- NA
      geno[i, , 1] <- a1; geno[i, , 2] <- a2
    }
    ind <- data.frame(id = ids, pop = paste0("pop", popn), time = NA_real_, use_for_ne = TRUE, stringsAsFactors = FALSE)
    ind$sample <- ind$pop
  } else if (fmt == "vcf") {
    eco_require("vcfR")
    v <- vcfR::read.vcfR(f, verbose = FALSE)
    gt <- vcfR::extract.gt(v, element = "GT")
    loci <- rownames(gt); ids <- colnames(gt)
    gt <- gsub("|", "/", gt, fixed = TRUE)
    a1 <- suppressWarnings(matrix(as.integer(sub("/.*", "", gt)), nrow(gt)))
    a2 <- suppressWarnings(matrix(as.integer(sub(".*/", "", gt)), nrow(gt)))
    geno <- array(NA_integer_, c(length(ids), length(loci), 2))
    geno[, , 1] <- t(a1) + 1L; geno[, , 2] <- t(a2) + 1L
    ind <- data.frame(id = ids, pop = "pop1", time = NA_real_, use_for_ne = TRUE, stringsAsFactors = FALSE)
    ind$sample <- ind$pop
  } else stop("unknown data format '", fmt, "' (use two_column, one_column, dosage, genepop or vcf)")
  if (!is.null(dd$popmap)) {
    pm <- pg_read_table(pg_resolve(config, dd$popmap))
    names(pm) <- tolower(names(pm))
    if (!all(c("id", "pop") %in% names(pm))) stop("popmap needs columns id and pop (optional: time, use_for_ne, sample)")
    k <- match(ind$id, pm$id)
    if (anyNA(k)) stop("popmap: ", sum(is.na(k)), " individuals in the genotype file are missing from the popmap, e.g. ", paste(utils::head(ind$id[is.na(k)], 3), collapse = ", "))
    ind$pop <- pm$pop[k]
    if (!is.null(pm$time)) ind$time <- as.numeric(pm$time[k])
    if (!is.null(pm$use_for_ne)) ind$use_for_ne <- toupper(pm$use_for_ne[k]) %in% c("TRUE", "T", "1", "YES", "Y")
    ind$sample <- if (!is.null(pm$sample)) pm$sample[k] else
      if (all(is.na(ind$time))) ind$pop else ifelse(is.na(ind$time), ind$pop, paste0(ind$pop, "_", ind$time))
  }
  dimnames(geno) <- list(ind$id, loci, NULL)
  kin <- NULL
  if (!is.null(dd$known_kin)) {
    kin <- pg_read_table(pg_resolve(config, dd$known_kin))
    names(kin) <- tolower(names(kin))
    if (!all(c("id1", "id2", "relationship") %in% names(kin))) stop("known_kin file needs columns id1, id2, relationship")
    kin$relationship <- tolower(gsub("[ -]", "", kin$relationship))
  }
  cn <- NULL
  if (!is.null(dd$census_N)) cn <- unlist(dd$census_N)
  marker <- if (!is.null(dd$marker)) dd$marker else {
    nal <- apply(geno, 2, function(g) length(unique(stats::na.omit(c(g)))))
    if (stats::median(nal) <= 2) "snp" else "msat"
  }
  list(geno = geno, ind = ind, known_kin = kin, census_N = cn,
       generation_time = if (is.null(dd$generation_time)) 1 else as.numeric(dd$generation_time),
       marker = marker, species = dd$species)
}

check_genotypes <- function(data) {
  pr <- character()
  g <- data$geno; ind <- data$ind
  if (length(dim(g)) != 3 || dim(g)[3] != 2) return("genotypes must be diploid (two alleles per locus)")
  if (anyDuplicated(ind$id)) pr <- c(pr, sprintf("individual IDs are duplicated: %s", paste(utils::head(unique(ind$id[duplicated(ind$id)]), 5), collapse = ", ")))
  half <- xor(is.na(g[, , 1]), is.na(g[, , 2]))
  if (any(half)) pr <- c(pr, sprintf("%d genotypes have only one allele scored (e.g. individual %s, locus %s); score both alleles or code both as missing",
                                     sum(half), dimnames(g)[[1]][which(half, arr.ind = TRUE)[1, 1]], dimnames(g)[[2]][which(half, arr.ind = TRUE)[1, 2]]))
  ms <- rowMeans(is.na(g[, , 1]))
  if (any(ms > 0.5)) pr <- c(pr, sprintf("%d individuals are missing more than half of their loci (e.g. %s); remove them before estimating Ne",
                                         sum(ms > 0.5), paste(utils::head(ind$id[ms > 0.5], 3), collapse = ", ")))
  ml <- colMeans(is.na(g[, , 1]))
  if (any(ml > 0.5)) pr <- c(pr, sprintf("%d loci are missing in more than half of individuals (e.g. %s); drop them",
                                         sum(ml > 0.5), paste(utils::head(dimnames(g)[[2]][ml > 0.5], 3), collapse = ", ")))
  ns <- table(ind$sample[ind$use_for_ne])
  if (any(ns < 10)) {
    small <- names(ns)[ns < 10]
    det <- vapply(small, function(s) {
      ids <- ind$id[ind$sample == s]
      near <- setdiff(names(ns)[ns >= 10][utils::adist(s, names(ns)[ns >= 10]) <= 3], s)
      sprintf("%s (%d: %s%s)%s", s, length(ids), paste(utils::head(ids, 3), collapse = ", "), if (length(ids) > 3) ", ..." else "",
              if (length(near)) sprintf(" - possibly a misspelling of %s?", near[1]) else "")
    }, "")
    pr <- c(pr, sprintf("sample(s) with fewer than 10 individuals: %s. Ne cannot be estimated from so few; if this is a typo in the site or year columns, correct it, otherwise these individuals are left out of the Ne estimates.",
                        paste(det, collapse = "; ")))
  }
  if (any(!is.na(ind$time)) && any(is.na(ind$time))) {
    w <- which(is.na(ind$time))
    pr <- c(pr, sprintf("%d individual(s) have no sampling time (e.g. %s); fill in the time column, otherwise they form a separate sample without a time and are not used by the temporal method",
                        length(w), paste(utils::head(ind$id[w], 3), collapse = ", ")))
  }
  if (!is.null(data$known_kin)) {
    bad <- setdiff(c(data$known_kin$id1, data$known_kin$id2), ind$id)
    if (length(bad)) pr <- c(pr, sprintf("known_kin lists IDs not in the genotype file: %s", paste(utils::head(bad, 5), collapse = ", ")))
  }
  pr
}
