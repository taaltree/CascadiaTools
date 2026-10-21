#!/usr/bin/env python3
"""Generate catalog.md (the module index agents choose from) from the
frontmatter of every modules/<id>/README.md. Rerun after adding a module.

Usage: python3 tools/make_catalog.py [--root ROOT]
"""
import argparse
import json
import re
from pathlib import Path

# Display order and grouping (matches the manuscript's Table 2).
DOMAINS = [
    ("Regression and causal inference", ["glm", "causal", "sem"]),
    ("Occupancy", ["occu-single", "occu-rota", "occu-waddle", "occu-spatial", "abund-interactions"]),
    ("Abundance", ["nmixture", "distance"]),
    ("Spatial capture-recapture", ["scr", "scr-smr", "scr-search", "scr-open"]),
    ("Demography and population dynamics", ["cjs", "multistate", "jolly-seber", "ipm", "timeseries", "pva"]),
    ("Distribution and community", ["sdm-presence-only", "jsdm", "community"]),
    ("Movement", ["movement-hmm", "movement-ssf"]),
    ("Conservation genetics", ["popgen-structure", "popgen-ne"]),
    ("Fisheries", ["fisheries-growth", "fisheries-stock"]),
]


def frontmatter(readme):
    text = readme.read_text()
    m = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    out = {}
    if m:
        for line in m.group(1).splitlines():
            if ":" in line:
                k, v = line.split(":", 1)
                out[k.strip()] = v.strip()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    root = Path(ap.parse_args().root)
    seen, lines = set(), [
        "# CascadiaTools module catalog",
        "",
        "Pick the module whose `when_to_use` matches the researcher's data and question, check its",
        "`not_for`, then read `modules/<id>/README.md` and `modules/<id>/VALIDATION.md` (known limits)",
        "before planning. Every module has a worked example in `modules/<id>/example/WALKTHROUGH.md`.",
        "If nothing fits, say so and follow the eco-new-module skill rather than forcing a poor fit.",
        "",
    ]
    index = []
    for dom, ids in DOMAINS:
        lines += [f"## {dom}", ""]
        for mid in ids:
            rd = root / "modules" / mid / "README.md"
            if not rd.is_file():
                continue
            fm = frontmatter(rd)
            seen.add(mid)
            lines += [f"### `{mid}`: {fm.get('title', mid)}", ""]
            for k, lab in [("when_to_use", "Use when"), ("not_for", "Not for"),
                           ("data_needed", "Data"), ("engine", "Engine")]:
                if fm.get(k):
                    lines.append(f"- **{lab}:** {fm[k]}")
            lines.append("")
            index.append({"id": mid, "domain": dom, **fm})
    extra = sorted(p.name for p in (root / "modules").iterdir()
                   if p.is_dir() and not p.name.startswith("_") and p.name not in seen)
    if extra:
        lines += ["## Other modules", ""] + [f"- `{e}`" for e in extra] + [""]
    (root / "catalog.md").write_text("\n".join(lines))
    (root / "docs" / "catalog.json").write_text(json.dumps(index, indent=1))
    print(f"catalog.md: {len(seen)} modules" + (f"; unlisted: {extra}" if extra else ""))


if __name__ == "__main__":
    main()
