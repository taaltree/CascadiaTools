---
name: cascadia-setup
description: Set up CascadiaTools on this computer (check R 4.5+, a C++ compiler for Nimble, and install the R packages; optionally Python for the GPU engines). Use the first time someone uses CascadiaTools, when a module says packages are missing, or when the user asks to install or update.
allowed-tools: Bash(Rscript *) Bash("${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" *)
---

# cascadia-setup — get a computer ready for CascadiaTools

The toolkit is at `${CLAUDE_PLUGIN_ROOT}`. Explain each step in one plain sentence;
many users have never used a terminal. Never install anything without saying what it
is and roughly how long it takes.

1. **R.** Run `"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" -e 'cat(R.version.string)'`.
   - Not found, or older than 4.5: tell the user to install the latest R from
     https://cran.r-project.org (Download R for Windows / macOS), with the default
     options, then restart Claude. Stop here until R works.
2. **Compiler** (needed for Bayesian Nimble models; most other modules work without it).
   - Windows: Rtools from https://cran.r-project.org/bin/windows/Rtools/ (match the R version).
   - macOS: run `xcode-select --install` and accept the dialog.
   - Linux: the distribution's build tools (`build-essential` or `gcc-c++`).
3. **Packages.** Run `"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/tools/install_packages.R"`
   in the background (20-40 minutes the first time). Tell the user it is safe to keep working.
   When it finishes, run it again with `--check` and report anything still missing,
   and which modules that affects.
4. **Test.** Run the template self-test:
   `"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/core/run.R" selftest glm --nsim 4`
   in a scratch folder in the user's project, and confirm it ends with a report path.
5. **Optional, GPU spatial capture-recapture:** only if the user has an NVIDIA GPU and
   large SCR problems: create a Python environment with the packages in
   `${CLAUDE_PLUGIN_ROOT}/tools/requirements-numpyro.txt` and set `CASCADIATOOLS_PYTHON`
   to its python. Everyone else can skip this; Nimble and secr cover the same models.
6. **Privacy.** Ask whether their data are sensitive (locations of rare species,
   data under agreements). If so, offer the `cascadia-private` skill.

Finish with one short paragraph: what is installed, what is not, and how to start
("Describe your data and your question, for example ...").
