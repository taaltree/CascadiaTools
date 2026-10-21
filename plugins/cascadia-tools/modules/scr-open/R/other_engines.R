# scr-open: NumPyro engine (NUTS, subprocess JSON in/out) and openCR
# (maximum likelihood, Efford & Schofield 2020) cross-check.

so_engine_numpyro <- function(prep, settings) {
  eco_require("jsonlite")
  py <- eco_python()
  if (!file.exists(py)) stop("NumPyro engine needs the Python env at ", py)
  td <- tempfile("scro_np_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ini <- so_inits_top(prep, prep$M)
  ent <- 1 - cumprod(1 - ini$gam); bi <- diff(c(0, ent / ent[length(ent)]))
  init <- list(p0 = ini$p0, sigma = ini$sigma, phi = ini$phi, psi = ent[length(ent)], b = I(as.numeric(bi)))
  if (!is.null(ini$bD)) init$bD <- I(as.numeric(ini$bD))
  inp <- list(d2 = prep$d2, Kjt = prep$usage * 1, XD = prep$XD, PD = ncol(prep$XD),
              nzi = I(prep$nzi - 1L), nzj = I(prep$nzj - 1L), nzt = I(prep$nzt - 1L), nzc = I(as.numeric(prep$nzc)),
              seen = prep$seen * 1L, n = prep$n, M = prep$M, sigma_max = prep$sigma_max,
              chains = settings$nchains, warmup = settings$numpyro_warmup, samples = settings$numpyro_samples,
              seed = settings$seed, x64 = !isTRUE(settings$numpyro_float32), init = init)
  fi <- file.path(td, "in.json"); fo <- file.path(td, "out.json")
  jsonlite::write_json(inp, fi, digits = NA, auto_unbox = TRUE, matrix = "rowmajor")
  env <- c("XLA_PYTHON_CLIENT_PREALLOCATE=false",
           if (!isTRUE(settings$numpyro_gpu)) "JAX_PLATFORMS=cpu",
           "XLA_FLAGS=--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
  t0 <- proc.time()[["elapsed"]]
  msg <- system2(py, c(shQuote(file.path(MODULE_DIR, "python", "scro_numpyro.py")), shQuote(fi), shQuote(fo)),
                 env = env, stdout = TRUE, stderr = TRUE)
  if (!file.exists(fo)) stop("NumPyro run failed:\n", paste(utils::tail(msg, 20), collapse = "\n"))
  res <- jsonlite::read_json(fo, simplifyVector = TRUE)
  dr <- res$draws; C <- settings$nchains
  chains <- lapply(seq_len(C), function(k) {
    m <- sapply(names(dr), function(v) as.numeric(matrix(dr[[v]], nrow = C)[k, ]))
    if (!is.matrix(m)) m <- matrix(m, nrow = 1, dimnames = list(NULL, names(dr)))
    m
  })
  list(chains = chains, t_total = proc.time()[["elapsed"]] - t0, engine = "numpyro marginal", timing = res$timing)
}

# ---------------------------------------------------------------- openCR
# Spatial JSSA (type JSSAsecrf: per-capita recruitment f_t, survival phi,
# static activity centres) fitted by maximum likelihood on the same mask.
# openCR's detection functions are hazard-based (HHN: p = 1 - exp(-lambda0
# exp(-d^2 / 2 sigma^2))), so p0 and sigma are not on exactly the same scale
# as the half-normal used here; N_t, phi and f_t are directly comparable.
so_engine_opencr <- function(prep, settings) {
  eco_require(c("secr", "openCR"))
  suppressPackageStartupMessages({ requireNamespace("secr"); requireNamespace("openCR") })
  f <- 1000
  K <- prep$K; TT <- prep$TT
  tdf <- data.frame(x = prep$txy[, 1] * f, y = prep$txy[, 2] * f); rownames(tdf) <- prep$traps$id
  tr <- secr::read.traps(data = tdf, detector = "proximity")
  # usage J x (K*T): first K_jt occasions of each period operating
  us <- matrix(0L, prep$J, K * TT)
  for (t in seq_len(TT)) for (j in seq_len(prep$J)) if (prep$usage[j, t] > 0) us[j, (t - 1) * K + seq_len(prep$usage[j, t])] <- 1L
  secr::usage(tr) <- us
  cap <- prep$captures
  ch <- secr::make.capthist(data.frame(Session = 1, ID = cap$individual, Occasion = (cap$period - 1) * K + cap$occasion,
                                       TrapID = cap$trap), tr, fmt = "trapID", noccasions = K * TT)
  attr(ch, "intervals") <- rep(c(rep(0, K - 1), 1), TT)[seq_len(K * TT - 1)]
  msk <- secr::read.mask(data = data.frame(x = prep$mxy[, 1] * f, y = prep$mxy[, 2] * f), spacing = prep$pix * f)
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(openCR::openCR.fit(ch, type = "JSSAsecrf", mask = msk, detectfn = "HHN",
                                             model = list(lambda0 ~ 1, sigma ~ 1, phi ~ 1, f ~ session),
                                             ncores = 1, trace = FALSE))
  pr <- stats::predict(fit)
  dv <- tryCatch(secr::derived(fit), error = function(e) { message("openCR derived() failed: ", conditionMessage(e)); NULL })
  list(fit = fit, predict = pr, derived = dv, t_total = proc.time()[["elapsed"]] - t0, engine = "openCR JSSAsecrf (ML, homogeneous density)")
}

# openCR results -> rows comparable with the module's estimates (point
# estimates for N_t; Wald intervals for phi, f_t, sigma).
so_opencr_rows <- function(prep, oc) {
  pr <- oc$predict; nm <- so_names(prep)
  pick <- function(p) { x <- pr[[p]]; if (is.data.frame(x)) x else do.call(rbind, x) }
  rows <- list()
  ph <- pick("phi")[1, ]
  rows[[1]] <- data.frame(param = "survival", estimate = ph$estimate, sd = ph$SE.estimate, lower = ph$lcl, upper = ph$ucl)
  sg <- pick("sigma")[1, ]
  rows[[2]] <- data.frame(param = "sigma", estimate = sg$estimate / 1000 * prep$cscale, sd = sg$SE.estimate / 1000 * prep$cscale,
                          lower = sg$lcl / 1000 * prep$cscale, upper = sg$ucl / 1000 * prep$cscale)
  ff <- pick("f")
  if (nrow(ff) >= prep$TT - 1) {
    ff <- ff[seq_len(prep$TT - 1), ]
    rows[[3]] <- data.frame(param = nm$f, estimate = ff$estimate, sd = ff$SE.estimate, lower = ff$lcl, upper = ff$ucl)
  }
  es <- oc$derived$estimates
  if (!is.null(es) && "D" %in% names(es)) {
    # D is animals per hectare of the mask (internal x1000 => metres)
    Nt <- es$D * prep$area * 100
    rows[[length(rows) + 1]] <- data.frame(param = nm$N, estimate = Nt[seq_len(prep$TT)], sd = NA, lower = NA, upper = NA)
  }
  do.call(rbind, rows)
}
