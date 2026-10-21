---
name: cascadia-start
description: Welcome and orientation for CascadiaTools - what analyses it can do, how to prepare data, how to try a worked example, and how a validated analysis works. Use when someone is new, asks what CascadiaTools is or can do, asks for help getting started, or wants to try an example.
---

# cascadia-start — orientation for new users

Toolkit location: `${CLAUDE_PLUGIN_ROOT}`. Module index: `${CLAUDE_PLUGIN_ROOT}/catalog.md`.

Speak plainly; assume the person is an ecologist, not a programmer.

1. **Check setup.** Run `"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/tools/install_packages.R" --check`.
   If anything is missing, offer the `cascadia-setup` skill first.
2. **Explain in four sentences** what happens in every analysis: you describe your data and
   question; the agent checks your files and writes a short plan; it simulates data like yours
   with known answers and confirms the method gets them right (planted truth); then it fits your
   data and gives you a report with plain-language results, publication figures and a methods
   paragraph. If the method cannot be trusted at your design, you are told so (OPEN), with why.
3. **What it covers:** 29 validated modules (count them in the catalog; never estimate), grouped as in the catalog: regression and causal inference; occupancy
   (single, multispecies, two-species, spatial, abundance-mediated interactions); abundance
   (N-mixture, distance sampling); spatial capture-recapture (cameras, hair snares, mark-resight,
   search tracks, open populations); survival and population dynamics (CJS, multistate,
   Jolly-Seber, integrated population models, trends, PVA); distributions and communities;
   movement; conservation genetics; fisheries. Ask what they work on and name the 1-3 modules
   that fit, with one sentence each from the catalog.
4. **Try an example.** Offer to run a module's worked example in a new folder
   `examples/<id>/` in their project: copy `${CLAUDE_PLUGIN_ROOT}/modules/<id>/example/` there,
   then follow its `WALKTHROUGH.md`, running the commands with the toolkit paths. This takes
   minutes for likelihood modules and longer for Bayesian ones; say which.
5. **Their own data.** Tell them what to bring (see the module README's "Data format"),
   suggest they keep originals in `data_raw/`, ask whether the data are sensitive
   (offer `cascadia-private`), then follow the `eco-analyze` skill.
6. **Study design.** If they have not collected data yet, the `eco-power` skill answers
   "how many sites / traps / animals / loci do I need?".

Citation for papers: Levi, T. CascadiaTools. Cascadia Wildlife Lab, Oregon State University.
https://github.com/taaltree/CascadiaTools (plus the packages each report lists).
