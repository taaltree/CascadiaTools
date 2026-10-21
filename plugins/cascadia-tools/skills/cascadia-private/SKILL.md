---
name: cascadia-private
description: Turn on data-private mode for this project so the AI agent works from data structure, check messages and results without reading raw records (sensitive species locations, data under agreements, human-subjects data). Use when the user says their data are sensitive, confidential or must not be shared, or asks how to keep data private.
---

# cascadia-private — analyse data without reading the raw records

**What it does.** The R toolkit runs on the user's computer. The AI model does not:
anything the agent reads or prints goes to the model provider. Private mode keeps raw
records out of the conversation. The agent sees only column names, dimensions, value
ranges, missing-value counts, `check_data()` messages and aggregated results.

## Turn it on

1. Ask which folders hold sensitive data. Default: `data_raw/` and `data/` in the project.
2. Write `.cascadia/private.json` in the project root:
   `{"private": true, "data_dirs": ["data_raw", "data"]}` (the folders the user named).
   The CascadiaTools guard hook reads this file and, while it exists, refuses tool calls
   that would show raw records: reading or searching those folders, and shell commands that
   print them. CascadiaTools commands (`core/run.R ...`) are allowed, because they print
   summaries and results, not records.
3. Add read-deny rules to the project's `.claude/settings.json` (create it if missing,
   merge if present), one per folder: `"permissions": {"deny": ["Read(./data_raw/**)", "Read(./data/**)"]}`.
   These block the Read, Grep and Glob tools and common file commands for those paths.
4. Tell the user, plainly, what this does and does not protect:
   - Protected: the agent cannot open or print the files; it works through summaries.
   - Still sent: summaries, column names, error messages, estimates and figures.
   - Not a guarantee: these are pattern-based guards, not a security boundary. A
     command written to evade them could still print data. For the strictest cases,
     run the toolkit without an AI model (`Rscript core/run.R ...` and the protocols
     as a written checklist), or use coordinates shifted and rotated before analysis.
   - Also check the model provider's data terms for their account or institution.

## While private mode is on

- Describe data with summaries only, e.g.
  `"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" -e 'd <- read.csv("data_raw/x.csv"); str(d[0, ]); print(dim(d)); print(colSums(is.na(d)))'`
  is still blocked (it names a private folder): use the module's `check_data` through
  `core/run.R fit` instead, which reports problems without printing records.
- Never paste rows, coordinates or identifiers into the conversation, a report or a commit.
- Reports go to the project's `results/` folder; remind the user to review them before sharing.

## Turn it off

Delete `.cascadia/private.json` and the deny rules, only when the user asks.
