# CascadiaTools

**Validated ecological analysis for AI coding agents.** CascadiaTools lets an AI agent
(Claude Code) run wildlife, fisheries, occupancy, abundance, movement, causal and
conservation-genetic analyses for you, using statistical code that has been written once,
checked by simulation and independently verified, instead of new, unchecked code every time.

You describe your data and your question in plain language. The agent:

1. checks your files and describes them back to you;
2. writes a short analysis plan;
3. simulates data like yours with known answers and confirms the method recovers them
   (*planted truth*); if it can't, you are told, with numbers;
4. fits your data, checks convergence and goodness of fit;
5. gives you a report with plain-language results, publication-ready figures and a
   methods paragraph.

Made by [Taal Levi](https://levi.cascadiawildlifelab.org) at the
[Cascadia Wildlife Lab](https://cascadiawildlifelab.org), Oregon State University.
Beginner's guide, examples and validation results: **https://cascadiawildlifelab.org/tools/**

## What it covers (29 modules)

| Area | Modules |
|---|---|
| Regression and causal inference | GLMs/GLMMs with biological interpretation; causal effects from DAGs; structural equation models |
| Occupancy | single-season and dynamic; multispecies (Rota); two-species conditional (Waddle); spatial (spOccupancy); abundance-mediated interactions |
| Abundance | N-mixture; distance sampling |
| Spatial capture-recapture | closed SCR (camera traps, hair snares, Poisson/binomial, masks, covariates); mark-resight and generalized mark-resight; search-encounter (sign search along tracks); open-population SCR |
| Demography | Cormack-Jolly-Seber; multistate; Jolly-Seber; integrated population models; population trends; population viability analysis |
| Distribution and community | presence-only SDMs; joint species distribution models (HMSC); ordination, PERMANOVA, beta diversity |
| Movement | hidden Markov models; step-selection analysis and home ranges |
| Conservation genetics | diversity and population structure; effective population size and relatedness |
| Fisheries | growth (von Bertalanffy) and condition; catch curves, stock-recruitment, surplus production |

Every module has a worked example (`plugins/cascadia-tools/modules/<module>/example/WALKTHROUGH.md`)
and a validation record (`VALIDATION.md`) with its recovery results and known limits.

## Install (about 30 minutes, once)

1. **Install R 4.5 or newer** from https://cran.r-project.org.
   Windows: also install Rtools. Mac: run `xcode-select --install` in Terminal.
2. **Install Claude Code**: the Claude desktop app (Code tab) or the terminal version;
   see https://claude.com/claude-code.
3. **Add CascadiaTools.** In Claude Code, type:
   ```
   /plugin marketplace add taaltree/CascadiaTools
   /plugin install cascadia-tools@cascadia
   ```
4. **Finish setup.** Open a folder for your project in Claude Code and type
   *"Set up CascadiaTools"*. It installs the R packages and runs a test.

Then try: *"Show me what CascadiaTools can do"*, or *"Run the occupancy example"*,
or describe your own data and question.

**Without the plugin system:** clone this repository and in Claude Code type
`/plugin marketplace add /path/to/CascadiaTools`, then install as above.
**Without any AI:** `Rscript plugins/cascadia-tools/core/run.R list` shows the modules;
each module's walkthrough lists the commands.

## Your data and privacy

The statistical code runs on your computer. The AI model does not: whatever the agent
reads or prints is sent to the model provider. For sensitive data, type
*"Turn on private mode"*: the agent then works from summaries and results without reading
raw records (see the `cascadia-private` skill for exactly what this does and does not protect).
For the most sensitive data, use the toolkit without AI.

## Citing

Please cite CascadiaTools (see `CITATION.cff`) and the packages each report lists.

## Licence

Code: MIT (`LICENSE`). Documentation, protocols and examples: CC BY 4.0 (`LICENSE-docs`).
