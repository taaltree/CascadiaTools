#!/usr/bin/env python3
"""Build Appendix S1 (worked examples for every module) from the per-module
summaries written by each module's worked example:

    results/<id>/si/si_summary.json   (+ the figure PNGs it points to)

Output: docs/si/appendix_S1.html with figures copied to docs/si/figs/<id>/.
The page is self-describing and regenerates whenever an example is rerun.

Usage: python3 tools/build_si.py [--root ROOT]
"""
import argparse
import html
import json
import shutil
from pathlib import Path

# Module order and grouping follow manuscript Table 2.
DOMAINS = [
    ("Regression and causal inference", ["glm", "causal", "sem"]),
    ("Occupancy", ["occu-single", "occu-rota", "occu-waddle", "occu-spatial", "abund-interactions"]),
    ("Abundance", ["nmixture", "distance"]),
    ("Spatial capture–recapture", ["scr", "scr-smr", "scr-search", "scr-open"]),
    ("Demography and population dynamics", ["cjs", "multistate", "jolly-seber", "ipm", "timeseries", "pva"]),
    ("Distribution and community", ["sdm-presence-only", "jsdm", "community"]),
    ("Movement", ["movement-hmm", "movement-ssf"]),
    ("Conservation genetics", ["popgen-structure", "popgen-ne"]),
    ("Fisheries", ["fisheries-growth", "fisheries-stock"]),
]
CORE = {"id", "title", "question", "prompt", "design", "truth", "recovery", "estimates",
        "diagnostics", "interpretation", "figures", "key_figure", "verdict", "limits", "commands"}

CSS = """
:root{--bg:#fbfaf7;--fg:#1d2321;--muted:#5d6762;--accent:#2f6f4f;--line:#d9ddd6;--card:#fff}
body{background:var(--bg);color:var(--fg);font:15px/1.55 'Source Sans 3',system-ui,sans-serif;margin:0}
main{max-width:1000px;margin:0 auto;padding:24px 16px 80px}
h1{font-size:1.8rem;margin:.2em 0}h2{font-size:1.45rem;margin-top:2.4em;border-bottom:2px solid var(--accent);padding-bottom:.2em}
h3{font-size:1.15rem;margin-top:1.8em}h4{font-size:1rem;margin:1.2em 0 .3em;color:var(--muted)}
.meta{color:var(--muted)}.tw{overflow-x:auto}
table{border-collapse:collapse;font-size:.86rem;font-variant-numeric:tabular-nums;margin:.4em 0}
th,td{border-bottom:1px solid var(--line);padding:3px 9px;text-align:right;vertical-align:top}
th:first-child,td:first-child{text-align:left}th{color:var(--muted);font-weight:600}
blockquote{margin:.6em 0;padding:8px 14px;background:var(--card);border-left:4px solid var(--accent);border-radius:4px}
figure{margin:1em 0}figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:4px;background:#fff}
figcaption{font-size:.88rem;color:var(--muted)}
.pass{color:#2f8f4f;font-weight:600}.fail{color:#c0392b;font-weight:600}
.verdict{display:inline-block;padding:2px 10px;border-radius:12px;font-weight:600}
.DONE{background:#e3f1e8;color:#1f6b3d}.OPEN{background:#fbeaea;color:#9b2c2c}
nav ol{columns:2;font-size:.92rem}code,pre{font-size:.85em}pre{background:var(--card);border:1px solid var(--line);padding:8px;overflow-x:auto}
"""


def web_copy(src, dst, width=1200):
    """Copy a figure for the web: print PNGs are 600 dpi, so downscale to
    `width` px (keeps aspect). Full-resolution files stay in results/."""
    try:
        from PIL import Image
        with Image.open(src) as im:
            if im.width > width:
                im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
            im.convert("RGB").save(dst.with_suffix(".webp"), "WEBP", quality=88, method=6)
        return dst.with_suffix(".webp").name
    except Exception:
        shutil.copy2(src, dst)
        return dst.name


def esc(x):
    return html.escape(str(x))


def fmt(v):
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.4g}"
    if v is None:
        return "–"
    if isinstance(v, dict):
        return "; ".join(f"{esc(k)} = {fmt(x)}" for k, x in v.items())
    if isinstance(v, list):
        return ", ".join(fmt(x) for x in v)
    return esc(v)


def kv_table(d):
    rows = "".join(f"<tr><td>{esc(k)}</td><td>{fmt(v)}</td></tr>" for k, v in d.items())
    return f"<div class='tw'><table><tr><th>Item</th><th>Value</th></tr>{rows}</table></div>"


def rows_table(rows, truth=None):
    if not rows:
        return "<p class='meta'>(none)</p>"
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    # Truth column only when most parameter names match the truth's names;
    # otherwise the truth table shown above the estimates is the reference.
    add_truth = truth is not None and "truth" not in cols and "param" in cols and \
        sum(r.get("param") in truth for r in rows) >= 0.5 * len(rows)
    head = "".join(f"<th>{esc(c)}</th>" for c in cols) + ("<th>truth</th>" if add_truth else "")
    body = []
    for r in rows:
        cells = []
        for c in cols:
            v = r.get(c)
            if c == "pass" and isinstance(v, bool):
                cells.append(f"<td class='{'pass' if v else 'fail'}'>{'PASS' if v else 'FAIL'}</td>")
            else:
                cells.append(f"<td>{fmt(v)}</td>")
        if add_truth:
            cells.append(f"<td>{fmt(truth.get(r.get('param')))}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<div class='tw'><table><tr>{head}</tr>{''.join(body)}</table></div>"


def anything(v, truth=None):
    if isinstance(v, dict):
        return kv_table(v)
    if isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
        return rows_table(v, truth)
    if isinstance(v, list):
        return "<ul>" + "".join(f"<li>{esc(x)}</li>" for x in v) + "</ul>"
    return f"<p>{esc(v)}</p>"


def section(mod_id, d, root, out_fig_dir, num):
    si = root / "results" / mod_id / "si"
    title = d.get("title", mod_id)
    verdict = str(d.get("verdict", "OPEN")).upper()
    verdict = "DONE" if verdict.startswith("DONE") else "OPEN"
    reason = d.get("verdict_reason") or d.get("verdict_note") or d.get("verdict_reasons")
    parts = [f"<h3 id='{esc(mod_id)}'>S1.{num} <code>{esc(mod_id)}</code>: {esc(title)}</h3>",
             f"<p><span class='verdict {verdict}'>{verdict}</span>"
             + (f" <span class='meta'>{esc(reason if isinstance(reason, str) else json.dumps(reason))}</span>" if reason else "")
             + "</p>",
             f"<h4>Question</h4><p>{esc(d.get('question', ''))}</p>"]
    if d.get("prompt"):
        parts.append(f"<h4>What the researcher typed</h4><blockquote>{esc(d['prompt'])}</blockquote>")
    if d.get("design"):
        parts.append("<h4>Study design</h4>" + anything(d["design"]))
    truth = d.get("truth") if isinstance(d.get("truth"), dict) else None
    if truth:
        parts.append("<h4>True values used to simulate the example data</h4>" + kv_table(truth))
    if d.get("recovery"):
        parts.append("<h4>Planted-truth recovery at this design</h4>" + anything(d["recovery"]))
    if d.get("estimates"):
        parts.append("<h4>Estimates from the example data</h4>" + anything(d["estimates"], truth))
    if d.get("diagnostics"):
        parts.append("<h4>Convergence and goodness of fit</h4>" + anything(d["diagnostics"]))
    if d.get("interpretation"):
        interp = d["interpretation"]
        interp = interp if isinstance(interp, list) else [interp]
        parts.append("<h4>Interpretation from the report</h4>" +
                     "".join(f"<blockquote>{esc(s)}</blockquote>" for s in interp))
    figs = d.get("figures") or []
    key = d.get("key_figure")
    figs = sorted(figs, key=lambda f: 0 if isinstance(f, dict) and f.get("id") == key else 1)
    if figs:
        parts.append("<h4>Figures</h4>")
        dest = out_fig_dir / mod_id
        dest.mkdir(parents=True, exist_ok=True)
        for f in figs:
            if not isinstance(f, dict):
                continue
            src = Path(f.get("file", ""))
            cands = [src, si / src, root / src, si / "fit" / "figures" / src.name]
            src = next((c for c in cands if c.is_file()), None)
            cap = f.get("caption", "")
            if src is None:
                parts.append(f"<p class='meta'>Figure {esc(f.get('id'))}: file not found.</p>")
                continue
            name = web_copy(src, dest / src.name)
            parts.append(f"<figure><img loading='lazy' src='figs/{esc(mod_id)}/{esc(name)}' "
                         f"alt='{esc(f.get('id', ''))}'><figcaption><b>{esc(f.get('id', ''))}.</b> {esc(cap)}</figcaption></figure>")
    if d.get("limits"):
        parts.append("<h4>Limits that matter at this design</h4>" + anything(d["limits"]))
    extra = [k for k in d if k not in CORE and not k.startswith("verdict") and k not in ("truth_source",)]
    for k in extra:
        parts.append(f"<h4>{esc(k.replace('_', ' '))}</h4>" + anything(d[k], truth))
    if d.get("commands"):
        parts.append("<h4>Commands to reproduce</h4><pre>" + esc("\n".join(d["commands"])) + "</pre>")
    walk = root / "modules" / mod_id / "example" / "WALKTHROUGH.md"
    if walk.is_file():
        parts.append(f"<p class='meta'>Step-by-step walkthrough: <code>modules/{esc(mod_id)}/example/WALKTHROUGH.md</code>.</p>")
    return "\n".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    a = ap.parse_args()
    root = Path(a.root)
    out = root / "docs" / "si"
    figdir = out / "figs"
    out.mkdir(parents=True, exist_ok=True)
    toc, body, missing, n = [], [], [], 0
    for dom, ids in DOMAINS:
        body.append(f"<h2>{esc(dom)}</h2>")
        for mid in ids:
            f = root / "results" / mid / "si" / "si_summary.json"
            if not f.is_file():
                missing.append(mid)
                continue
            n += 1
            d = json.loads(f.read_text())
            toc.append(f"<li><a href='#{esc(mid)}'>{esc(mid)}</a>: {esc(d.get('title', ''))}</li>")
            body.append(section(mid, d, root, figdir, n))
    intro = (
        "<h1>Appendix S1. Worked examples for every CascadiaTools module</h1>"
        "<p class='meta'>Supporting Information for: CascadiaTools: a validated, AI-guided analysis harness for "
        "ecology, wildlife, fisheries and conservation genetics. T. Levi, Cascadia Wildlife Lab.</p>"
        "<p>Each example follows the same workflow an AI agent follows for a researcher. It gives the question, the "
        "prompt the researcher typed, the study design, and the true values used to simulate the example dataset. "
        "It then shows planted-truth recovery at that design (the acceptance gate), estimates from the example data "
        "set against the truth, convergence and goodness of fit, the report's interpretation, and the module's "
        "default publication figures. A verdict of DONE means every interpreted estimand passed the gate at this "
        "design. OPEN means at least one did not; the reason is stated, and the example shows how the harness reports "
        "it. All example data are simulated, and every result can be reproduced with the listed commands.</p>"
        + (f"<p class='meta'>Not yet available: {esc(', '.join(missing))}.</p>" if missing else "")
        + f"<nav><h4>Contents</h4><ol>{''.join(toc)}</ol></nav>"
    )
    page = ("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
            "content='width=device-width,initial-scale=1'><title>CascadiaTools Appendix S1</title>"
            f"<style>{CSS}</style></head><body><main>{intro}{''.join(body)}</main></body></html>")
    (out / "appendix_S1.html").write_text(page)
    print(f"wrote {out / 'appendix_S1.html'}: {n} examples; missing: {missing}")


if __name__ == "__main__":
    main()
