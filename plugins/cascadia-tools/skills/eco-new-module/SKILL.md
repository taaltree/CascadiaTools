---
name: eco-new-module
description: Build and validate a new CascadiaTools analysis module (new model family or variant). Use when no existing module fits a researcher's question or when extending the toolkit.
---

# eco-new-module

Follow `${CLAUDE_PLUGIN_ROOT}/docs/MODULE_SPEC.md` exactly; start by copying
`${CLAUDE_PLUGIN_ROOT}/modules/_template/`. Order of work:

1. **Read the primary literature** for the model (download and read the
   paper, not memory) and write the model equations into README.md before
   coding. Note each assumption.
2. Write `simulate()` straight from the equations: it is the planted truth.
   Keep it independent of the fitting code (no shared likelihood helpers).
3. Write `fit()`; summarise to contract estimates on biological scales.
4. `selftest`, then `recover` nsim >= 100 at the default design, and a harder
   design; negative controls; cross-check against an independent
   implementation (eco-acceptance-gate step 4).
5. `read_data`, `check_data`, `interpret`, `plots`, `methods_text`, and
   an `example/` dataset with `config.yaml`; run `fit` on it and read the
   HTML report.
6. Write VALIDATION.md (numbers, not adjectives; Known limits).
7. Add the module to `${CLAUDE_PLUGIN_ROOT}/catalog.md` (one line: id, question, data).
