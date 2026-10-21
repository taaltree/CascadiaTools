# jsdm: read researchers' community data and check it.
#
# Files (CSV):
#   community: species records. Either
#     long : one row per site x species record: site, species, value
#            (count or 0/1). Site x species combinations that are not listed
#            are zeros (absences) for every site in the sites file. Repeated
#            rows (several visits) are summed (counts) / combined (presence).
#     wide : one row per site, one column per species (site column + species columns)
#   sites: one row per sampling unit: site, covariates [, x, y coordinates]
#   traits (optional): one row per species: species, trait columns
#            [, taxonomy columns such as order, family, genus]
#   phylogeny (optional): Newick tree file whose tip labels are the species
#            names (spaces may be written as underscores)

# "did you mean" hints for misspelled ids (closest names by edit distance)
jsdm_suggest <- function(bad, pool, k = 3) {
  bad <- utils::head(bad, k)
  vapply(bad, function(b) {
    d <- utils::adist(b, pool, ignore.case = TRUE)[1, ]
    if (!length(d) || min(d) > max(2, nchar(b) %/% 3)) sprintf("'%s'", b)
    else sprintf("'%s' (did you mean '%s'?)", b, pool[which.min(d)])
  }, "")
}

read_jsdm <- function(config) {
  dd <- config$data; dir <- config$.dir %||% "."
  rd <- function(f) utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE, check.names = FALSE,
                                    strip.white = TRUE)
  sc <- dd$site_col %||% "site"; spc <- dd$species_col %||% "species"
  st <- rd(dd$sites)
  if (!sc %in% names(st))
    stop(sprintf("the sites file (%s) has no column '%s'. Its columns are: %s. Set 'data: site_col:' in the config.",
                 dd$sites, sc, paste(names(st), collapse = ", ")), call. = FALSE)
  st[[sc]] <- as.character(st[[sc]])
  com <- rd(dd$community)
  fmt <- dd$format %||% (if (spc %in% names(com)) "long" else "wide")
  if (!sc %in% names(com))
    stop(sprintf("the community file (%s) has no column '%s'. Its columns are: %s.",
                 dd$community, sc, paste(names(com), collapse = ", ")), call. = FALSE)
  com[[sc]] <- as.character(com[[sc]])
  notes <- character()
  if (fmt == "long") {
    vc <- dd$value_col %||% "count"
    for (k in c(spc, vc)) if (!k %in% names(com))
      stop(sprintf("the community file (%s) has no column '%s'. Its columns are: %s. Set 'data: species_col:' / 'value_col:' in the config.",
                   dd$community, k, paste(names(com), collapse = ", ")), call. = FALSE)
    v <- suppressWarnings(as.numeric(com[[vc]]))
    if (any(is.na(v) & !is.na(com[[vc]]) & nzchar(as.character(com[[vc]]))))
      stop(sprintf("column '%s' in %s has non-numeric values, e.g. '%s'. Use numbers (counts or 0/1).",
                   vc, dd$community, com[[vc]][which(is.na(v))[1]]), call. = FALSE)
    com$.v <- v
    unk <- setdiff(unique(com[[sc]]), st[[sc]])
    if (length(unk))
      stop(sprintf("%d site(s) in the community file are missing from the sites file: %s. Check the spelling of the site ids in both files; every surveyed site needs a row (with covariates) in the sites file.",
                   length(unk), paste(jsdm_suggest(unk, st[[sc]]), collapse = ", ")), call. = FALSE)
    if (anyDuplicated(com[, c(sc, spc)]))
      notes <- c(notes, "Some site x species combinations appear on several rows (e.g. repeated visits); they were summed per site.")
    spp <- sort(unique(as.character(com[[spc]])))
    Y <- matrix(0, nrow(st), length(spp), dimnames = list(st[[sc]], spp))
    agg <- stats::aggregate(com$.v, list(site = com[[sc]], species = as.character(com[[spc]])), sum, na.rm = TRUE)
    Y[cbind(match(agg$site, st[[sc]]), match(agg$species, spp))] <- agg$x
  } else {
    spcols <- dd$species_cols %||% setdiff(names(com), c(sc, dd$exclude_cols))
    unk <- setdiff(com[[sc]], st[[sc]])
    if (length(unk))
      stop(sprintf("%d site(s) in the community file are missing from the sites file: %s. Check the spelling of the site ids in both files.",
                   length(unk), paste(jsdm_suggest(unk, st[[sc]]), collapse = ", ")), call. = FALSE)
    Y <- as.matrix(com[match(st[[sc]], com[[sc]]), spcols, drop = FALSE])
    rownames(Y) <- st[[sc]]
    if (any(is.na(match(st[[sc]], com[[sc]]))))
      stop("some sites in the sites file have no row in the wide community file; add them (all zeros if nothing was found).", call. = FALSE)
  }
  resp <- dd$response %||% (if (all(Y %in% c(0, 1))) "presence" else "count")
  if (resp == "presence") Y <- (Y > 0) * 1
  # drop species never recorded (nothing to estimate) and, optionally, very rare ones
  minocc <- dd$min_occurrences %||% 1
  occ <- colSums(Y > 0)
  drop <- names(occ)[occ < minocc]
  if (length(drop)) {
    notes <- c(notes, sprintf("%d species recorded at fewer than %d site(s) were left out of the model: %s.",
                              length(drop), minocc, paste(drop, collapse = ", ")))
    Y <- Y[, occ >= minocc, drop = FALSE]
  }
  xf <- dd$formula %||% (if (!is.null(dd$covariates)) paste("~", paste(dd$covariates, collapse = " + ")) else "~ 1")
  tr <- NULL; trf <- NULL; tree <- NULL
  labs <- stats::setNames(colnames(Y), colnames(Y))
  if (!is.null(dd$traits)) {
    tr <- rd(dd$traits)
    tsc <- dd$trait_species_col %||% spc
    if (!tsc %in% names(tr))
      stop(sprintf("the traits file (%s) has no species column '%s'. Its columns are: %s.",
                   dd$traits, tsc, paste(names(tr), collapse = ", ")), call. = FALSE)
    tr[[tsc]] <- as.character(tr[[tsc]])
    miss <- setdiff(colnames(Y), tr[[tsc]])
    if (length(miss))
      stop(sprintf("%d species have no row in the traits file: %s. Check the spelling in both files, add the missing rows, or drop traits from the config.",
                   length(miss), paste(jsdm_suggest(miss, tr[[tsc]]), collapse = ", ")), call. = FALSE)
    tr <- tr[match(colnames(Y), tr[[tsc]]), , drop = FALSE]
    names(tr)[names(tr) == tsc] <- "species"
    trf <- dd$trait_formula %||% (if (!is.null(dd$trait_cols)) paste("~", paste(dd$trait_cols, collapse = " + ")) else NULL)
    if (!is.null(dd$label_col) && !dd$label_col %in% names(tr))
      stop(sprintf("label_col '%s' is not a column of the traits file (%s). Its columns are: %s.",
                   dd$label_col, dd$traits, paste(names(tr), collapse = ", ")), call. = FALSE)
    if (!is.null(dd$label_col) && dd$label_col %in% names(tr)) labs <- stats::setNames(tr[[dd$label_col]], colnames(Y))
    if (!is.null(dd$taxonomy)) {
      eco_require("ape")
      badtx <- setdiff(dd$taxonomy, names(tr))
      if (length(badtx))
        stop(sprintf("taxonomy column(s) %s are not columns of the traits file (%s). Its columns are: %s.",
                     paste0("'", badtx, "'", collapse = ", "), dd$traits, paste(names(tr), collapse = ", ")), call. = FALSE)
      tx <- tr[, c(dd$taxonomy, "species"), drop = FALSE]
      tx[] <- lapply(tx, function(v) factor(gsub("[^A-Za-z0-9]+", "_", v)))
      fml <- stats::as.formula(paste("~", paste(c(dd$taxonomy, "species"), collapse = "/")))
      tree <- ape::as.phylo(fml, data = tx, collapse = FALSE)
      tree$edge.length <- rep(1, nrow(tree$edge))
      # restore original species names on the tips
      tree$tip.label <- colnames(Y)[match(tree$tip.label, gsub("[^A-Za-z0-9]+", "_", colnames(Y)))]
      notes <- c(notes, sprintf("Taxonomic tree built from %s (all branch lengths 1), used as the phylogenetic correlation matrix.",
                                paste(dd$taxonomy, collapse = " / ")))
    }
  }
  if (!is.null(dd$phylogeny)) {
    eco_require("ape")
    tree <- ape::read.tree(file.path(dir, dd$phylogeny))
    map <- colnames(Y)[match(tree$tip.label, colnames(Y))]
    alt <- colnames(Y)[match(gsub("_", " ", tree$tip.label), colnames(Y))]
    map[is.na(map)] <- alt[is.na(map)]
    lost <- setdiff(colnames(Y), map)
    if (length(lost))
      stop(sprintf("%d species are not tips of the phylogeny, e.g. %s.", length(lost),
                   paste(utils::head(lost, 5), collapse = ", ")), call. = FALSE)
    lookup <- stats::setNames(map, tree$tip.label)
    tree <- ape::keep.tip(tree, tree$tip.label[!is.na(map)])
    tree$tip.label <- unname(lookup[tree$tip.label])
  }
  coords <- NULL
  if (!is.null(dd$x_col) && !is.null(dd$y_col)) {
    fac <- if (identical(tolower(dd$coord_units %||% "m"), "m")) 1000 else 1
    coords <- cbind(x = st[[dd$x_col]] / fac, y = st[[dd$y_col]] / fac)
  }
  list(Y = Y, site_data = st, trait_data = tr, x_formula = xf, tr_formula = trf,
       response = resp, coords = coords, tree = tree, species_labels = labs,
       read_notes = notes, labels = dd$labels %||% config$figures$labels)
}

check_jsdm <- function(data) {
  pr <- character()
  Y <- data$Y
  if (ncol(Y) < 2) pr <- c(pr, "fewer than 2 species: use a single-species model (glm or occu-single).")
  if (nrow(Y) < 20) pr <- c(pr, sprintf("only %d sites: a JSDM needs at least ~20 (better 50+) sampling units to estimate residual associations.", nrow(Y)))
  if (anyNA(Y)) pr <- c(pr, "the community matrix has missing values; enter 0 for species not found at surveyed sites.")
  if (data$response == "count" && any(Y < 0 | Y != round(Y), na.rm = TRUE))
    pr <- c(pr, "counts must be non-negative whole numbers (for cover/biomass data use a different model).")
  vars <- all.vars(stats::as.formula(data$x_formula))
  miss <- setdiff(vars, names(data$site_data))
  if (length(miss)) pr <- c(pr, sprintf("covariate(s) %s are not columns of the sites file. Its columns are: %s.",
                                        paste(jsdm_suggest(miss, names(data$site_data), k = length(miss)), collapse = ", "),
                                        paste(names(data$site_data), collapse = ", ")))
  for (v in intersect(vars, names(data$site_data))) {
    x <- data$site_data[[v]]
    if (anyNA(x)) pr <- c(pr, sprintf("covariate '%s' has %d missing value(s) (sites %s); fill or drop those sites.",
                                      v, sum(is.na(x)), paste(utils::head(data$site_data[[1]][is.na(x)], 5), collapse = ", ")))
    if (length(unique(x[!is.na(x)])) < 2) pr <- c(pr, sprintf("covariate '%s' does not vary across sites.", v))
  }
  if (!is.null(data$tr_formula)) {
    tv <- all.vars(stats::as.formula(data$tr_formula))
    miss <- setdiff(tv, names(data$trait_data))
    if (length(miss)) pr <- c(pr, sprintf("trait(s) %s are not columns of the traits file. Its columns are: %s.",
                                          paste(jsdm_suggest(miss, names(data$trait_data), k = length(miss)), collapse = ", "),
                                          paste(names(data$trait_data), collapse = ", ")))
    for (v in intersect(tv, names(data$trait_data)))
      if (anyNA(data$trait_data[[v]])) pr <- c(pr, sprintf("trait '%s' is missing for species %s.", v,
                                                            paste(data$trait_data$species[is.na(data$trait_data[[v]])], collapse = ", ")))
  }
  occ <- colSums(Y > 0)
  if (data$response == "presence") {
    ubiq <- names(occ)[occ == nrow(Y)]
    if (length(ubiq)) pr <- c(pr, sprintf("species present at every site carry no information on responses: %s. Remove them (or model counts).",
                                          paste(ubiq, collapse = ", ")))
  }
  pr
}
