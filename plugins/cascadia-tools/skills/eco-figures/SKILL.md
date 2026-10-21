---
name: eco-figures
description: Produce and customise publication-ready figures for any CascadiaTools analysis (journal sizes, labels, units, colours, formats, multi-panel layouts). Use when a researcher wants figures for a paper, talk or report, or asks to change how a figure looks.
---

# eco-figures — publication figures without code

Every module has a figure suite (see `${CLAUDE_PLUGIN_ROOT}/docs/FIGURES.md`). `${CLAUDE_PLUGIN_ROOT}/core/run.R fit`
exports it automatically to `<out>/figures/` with a `gallery.html`. You
customise figures by editing the `figures:` block of the project config and
redrawing from the saved fit. Never refit just to change a figure.

```
"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/core/run.R" figures <module> --fit <out>/<module>_fit.rds --config config.yaml
```

## Conversation loop

1. **Ask only what you need:**
   - which journal or venue (sets widths and font size);
   - which figures go in the paper (`include:`);
   - covariate names and units they want shown;
   - species, state or group names;
   - colour or greyscale.
   Defaults are good for everything else.
2. **Translate requests into config, not code.**
   - "Show forest as % not z-scores" → `units: {forest: {center, scale}}` with the values from data intake, plus `labels: {forest: "Forest cover (%)"}`.
   - "Make it two-column" → `per_figure: {<id>: {width: double}}`.
   - "Black and white for print" → `greyscale: true`.
   - "Combine the three detection plots" → a small composition script in the project folder using `compose_figures()`.
3. **Redraw, then look.** Open the PNG with the Read tool and check it like
   a reviewer would:
   - labels overlap or are clipped;
   - text is too small at the final width (<6 pt);
   - units are missing;
   - legend covers data;
   - colours don't separate in greyscale;
   - an axis implies precision the data lack;
   - the y range hides the data or exaggerates an effect;
   - extrapolation beyond the observed covariate range.
   Fix what you find before showing the researcher.
4. **Deliver:**
   - vector PDF for the journal, TIFF if they require it, PNG for slides;
   - the edited captions;
   - the config block, so they can reproduce or tweak the figures later.

## Rules

- Don't change what a figure means while restyling it: the CI level, what
  is held constant, back-transformations. If the researcher asks for
  something misleading (truncating the axis to exaggerate an effect,
  removing the uncertainty band, a dual axis), explain why it misleads and
  offer an honest alternative.
- Figures from simulated data (recovery, power) are labelled as simulations.
- Project-specific custom figures go in `projects/<name>/figures.R`, built
  from the `fig_*` constructors and saved with `save_figure()` so they
  share the style. Don't edit module code for one paper.
