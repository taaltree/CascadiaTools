# CascadiaTools plugin

Validated analysis modules and working protocols that let Claude Code run ecological,
wildlife, fisheries and conservation-genetic analyses correctly: every analysis checks that
the method recovers known truth at your study design before your data are interpreted, and
ends with a plain-language report, publication figures and a methods paragraph.

Install in Claude Code:

```
/plugin install cascadia-tools --marketplace taaltree/CascadiaTools
```

(or `/plugin marketplace add taaltree/CascadiaTools` then `/plugin install cascadia-tools@cascadia`).

Start with *"Show me what CascadiaTools can do"* (skill `cascadia-start`), or
*"Set up CascadiaTools"* (skill `cascadia-setup`) on a new computer. Requires R 4.5+.

Skills: `eco-analyze` (end-to-end workflow), `eco-planted-truth`, `eco-acceptance-gate`,
`eco-analysis-plan`, `eco-data-intake`, `eco-convergence`, `eco-model-checking`,
`eco-interpretation`, `eco-power`, `eco-figures`, `eco-new-module`, `cascadia-start`,
`cascadia-setup`, `cascadia-private`. Module index: `catalog.md`.

Guide and examples: https://cascadiawildlifelab.org/tools/ · Cascadia Wildlife Lab, Oregon State University.
