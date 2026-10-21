# Self-contained HTML reports (plots embedded as base64 PNG) so a
# researcher gets one file they can open, email or archive.
#
# analysis_report(): real-data fit -> methods, estimates, diagnostics,
#   plain-language interpretation, module plots (conditional effects, ...).
# recovery_report(): simulation-recovery / power results with the gate.

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

html_table <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return("<p><em>(none)</em></p>")
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(v) formatC(signif(v, digits), format = "fg", flag = "#", digits = digits))
  head <- paste0("<tr>", paste0("<th>", html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  body <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", html_escape(trimws(r)), "</td>", collapse = ""), "</tr>"))
  paste0("<div class='tw'><table>", head, paste(body, collapse = ""), "</table></div>")
}

plot_png64 <- function(p, width = 7, height = 4.5, res = 110) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = width, height = height, units = "in", res = res)
  ok <- tryCatch({ print(p); TRUE }, error = function(e) { message("plot failed: ", conditionMessage(e)); FALSE })
  grDevices::dev.off()
  if (!ok) return("")
  eco_require("base64enc")
  sprintf("<figure><img alt='plot' src='data:image/png;base64,%s'/></figure>",
          base64enc::base64encode(f))
}

.report_css <- "
:root{--bg:#fbfaf7;--fg:#1d2321;--muted:#5d6762;--accent:#2f6f4f;--line:#d9ddd6;--card:#fff}
@media (prefers-color-scheme: dark){:root{--bg:#141816;--fg:#e4e8e5;--muted:#9aa59f;--accent:#7cc39c;--line:#2c3430;--card:#1b201d}}
body{background:var(--bg);color:var(--fg);font:16px/1.55 system-ui,-apple-system,Segoe UI,sans-serif;margin:0}
main{max-width:980px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:1.7rem;margin:.2em 0}h2{font-size:1.25rem;margin-top:2em;border-bottom:1px solid var(--line);padding-bottom:.2em}
.meta{color:var(--muted);font-size:.9rem}.tw{overflow-x:auto}
table{border-collapse:collapse;font-size:.88rem;font-variant-numeric:tabular-nums;margin:.6em 0}
th,td{border-bottom:1px solid var(--line);padding:4px 10px;text-align:right}th:first-child,td:first-child{text-align:left}
th{color:var(--muted);font-weight:600}
figure{margin:1em 0}img{max-width:100%;height:auto;background:#fff;border-radius:6px}
.pass{color:#2f8f4f;font-weight:600}.fail{color:#c0392b;font-weight:600}
.callout{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--accent);padding:10px 14px;border-radius:6px}
code{font-size:.9em}"

html_page <- function(title, sections, subtitle = NULL) {
  paste0("<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<meta name='viewport' content='width=device-width,initial-scale=1'>",
         "<title>", html_escape(title), "</title><style>", .report_css, "</style></head><body><main>",
         "<h1>", html_escape(title), "</h1>",
         if (!is.null(subtitle)) paste0("<p class='meta'>", html_escape(subtitle), "</p>"),
         paste(sections, collapse = "\n"), cascadia_footer(), "</main></body></html>")
}

cascadia_footer <- function() {
  paste0("<footer style='margin-top:3em;padding-top:1em;border-top:1px solid var(--line);font-size:.85rem;color:var(--muted)'>",
         "Made with <a href='", CASCADIA_URL, "'>", CASCADIA_NAME, "</a> · Cascadia Wildlife Lab, Oregon State University · ",
         "Please cite: ", html_escape(CASCADIA_CITE), "</footer>")
}

section <- function(h, ...) paste0("<h2>", html_escape(h), "</h2>", paste0(..., collapse = "\n"))
paras <- function(x) paste0("<p>", html_escape(x), "</p>", collapse = "\n")

analysis_report <- function(module, fit, data, settings, file, config = NULL) {
  est <- fit$estimates
  dg <- fit$diagnostics
  conv <- if (isTRUE(dg$converged)) "<span class='pass'>converged</span>"
          else "<span class='fail'>NOT converged: do not interpret these estimates</span>"
  sec <- list(
    section("Question", paras(module$description)),
    section("Methods", paras(if (is.function(module$methods_text)) module$methods_text(fit, data, settings)
                             else paste("Engine:", module$engine))),
    section("Diagnostics", "<p class='callout'>Model ", conv,
            if (!is.null(dg$rhat_max)) sprintf("; max R-hat %.3f", dg$rhat_max),
            if (!is.null(dg$ess_min)) sprintf("; min bulk ESS %.0f", dg$ess_min), "</p>",
            if (length(dg$notes)) paras(dg$notes)),
    section("Estimates", html_table(est)),
    if (is.function(module$interpret)) section("What this means biologically", paras(module$interpret(fit, data))),
    if (is.function(module$plots)) {
      pl <- tryCatch(module$plots(fit, data), error = function(e) list())
      section("Figures", paste(vapply(names(pl), function(n)
        paste0("<h3>", html_escape(n), "</h3>", plot_png64(pl[[n]])), ""), collapse = "\n"))
    },
    if (length(module$references)) section("References", paste0("<ul>", paste0("<li>", html_escape(module$references), "</li>", collapse = ""), "</ul>")),
    section("Reproduce", "<pre><code>", html_escape(paste0(
      "Rscript core/run.R fit ", module$id, if (!is.null(config)) paste0(" --config ", config))), "</code></pre>")
  )
  writeLines(html_page(module$title, unlist(sec), paste(CASCADIA_NAME, "·", module$id, "·", format(Sys.time(), "%Y-%m-%d %H:%M"))), file)
  invisible(file)
}

recovery_plot <- function(reps) {
  eco_require("ggplot2")
  d <- reps[!is.na(reps$param), ]
  d$covers <- d$lower <= d$truth & d$upper >= d$truth
  ggplot2::ggplot(d, ggplot2::aes(x = estimate)) +
    ggplot2::geom_histogram(bins = 25, fill = "#7fa99b", colour = "white") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = truth), colour = "#c0392b", linewidth = 0.9) +
    ggplot2::facet_wrap(~param, scales = "free") +
    ggplot2::labs(x = "Estimate across simulated datasets (red line = true value)", y = "Replicates") +
    ggplot2::theme_minimal(base_size = 11)
}

power_plot <- function(pw) {
  eco_require("ggplot2")
  d <- pw[!is.na(pw$power), ]
  ggplot2::ggplot(d, ggplot2::aes(value, power, colour = param)) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = 2, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = d$vary[1], y = "Power (95% interval excludes no-effect)", colour = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}

recovery_report <- function(module, res, file) {
  s <- res$summary
  g <- if (nrow(s)) recovery_gate(s) else NULL
  if (!is.null(g)) s$gate <- ifelse(g$pass, "PASS", "FAIL")
  verdict <- if (!is.null(g) && all(g$pass)) "<p class='callout'><span class='pass'>PASS</span>: every parameter recovered without detectable bias and with near-nominal interval coverage.</p>"
             else "<p class='callout'><span class='fail'>FAIL or incomplete</span>: at least one parameter is biased, under-covered or failed to converge at this design. See the table.</p>"
  sec <- list(
    section("What was simulated", paras(module$description),
            "<h3>True parameters</h3>", kv_table(res$params),
            "<h3>Design</h3>", kv_table(res$design)),
    section("Verdict", verdict, html_table(s)),
    section("Sampling distributions", plot_png64(recovery_plot(res$reps), 8, 5)),
    section("How to read this", paras(c(
      "Bias: average of (estimate - truth). Relative bias divides by the true value.",
      "Coverage: share of replicates whose 95% interval contains the true value; it should be near 0.95.",
      "Power: share of replicates whose 95% interval excludes the no-effect value.",
      sprintf("Each replicate used seed %d*100000 + replicate number; rerun any one with module$simulate(params, design, seed).", res$seed))))
  )
  writeLines(html_page(paste(module$title, "- simulation recovery"), unlist(sec),
                       sprintf("%d replicates · %s · %s", res$nsim, CASCADIA_NAME, format(Sys.time(), "%Y-%m-%d"))), file)
  invisible(file)
}

kv_table <- function(x) html_table(data.frame(name = names(x),
  value = vapply(x, function(v) paste(format(v), collapse = ", "), ""), row.names = NULL))

power_report <- function(module, pw, file, params, design, nsim) {
  vary <- pw$vary[1]
  sec <- list(
    section("Question", paras(module$description)),
    section("Simulation settings", paras(sprintf(
      "Power was estimated from %d simulated datasets at each value of '%s'. All other settings were held at the values below.", nsim, vary)),
      "<h3>True parameters</h3>", kv_table(params), "<h3>Design</h3>", kv_table(design)),
    section("Power curve", plot_png64(power_plot(pw)),
            paras("Dashed line: the conventional 80% power target. Power is the share of simulated datasets whose 95% interval excludes the no-effect value."),
            html_table(pw[, intersect(c("value", "param", "truth", "power", "coverage", "rel_bias", "rmse", "n_converged"), names(pw))])),
    section("Caveat", paras("Power is only meaningful where bias is small and coverage is near 95%. Check those columns at each design value before choosing a sample size."))
  )
  writeLines(html_page(paste(module$title, "- power analysis"), unlist(sec),
                       paste(CASCADIA_NAME, "·", format(Sys.time(), "%Y-%m-%d"))), file)
  invisible(file)
}

# Conditional-effects plot shared by all modules: response (on its natural
# scale) against one covariate, other covariates held at their means, 95%
# band, rug of observed covariate values. `pred` is a data.frame with
# columns x, estimate, lower, upper (and optional `group` for several
# curves, e.g. species or states). Modules compute pred however suits the
# engine (posterior draws, delta method, predict()).
# Fixed (non-mapped) aesthetics for a layer, only when there is no grouping.
# In ggplot2 4.x a layer parameter set to NULL (fill = NULL) is not "unset":
# it overrides the colour/fill mapping and drops the grouping, so grouped
# curves came out black and ribbons zig-zagged between groups.
eco_fixed <- function(has_g, ...) if (has_g) list() else list(...)

cond_effect_plot <- function(pred, xlab, ylab, observed_x = NULL, title = NULL) {
  eco_require("ggplot2")
  has_g <- "group" %in% names(pred)
  aes_line <- if (has_g) ggplot2::aes(x, estimate, colour = group, fill = group) else ggplot2::aes(x, estimate)
  p <- ggplot2::ggplot(pred, aes_line) +
    do.call(ggplot2::geom_ribbon, c(list(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.22, colour = NA),
                                    eco_fixed(has_g, fill = "#2f6f4f"))) +
    do.call(ggplot2::geom_line, c(list(linewidth = 1), eco_fixed(has_g, colour = "#2f6f4f"))) +
    ggplot2::labs(x = xlab, y = ylab, title = title, colour = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 12)
  if (!is.null(observed_x))
    p <- p + ggplot2::geom_rug(data = data.frame(x = observed_x), ggplot2::aes(x = x),
                               inherit.aes = FALSE, alpha = 0.4, sides = "b")
  p
}
