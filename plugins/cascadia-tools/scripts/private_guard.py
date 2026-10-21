#!/usr/bin/env python3
"""CascadiaTools data-private guard (PreToolUse hook).

Does nothing unless the project has .cascadia/private.json with "private": true.
Then it denies tool calls that would show raw records from the listed data folders:
Read/Grep/Glob on those paths, and Bash commands that mention them, except
CascadiaTools pipeline commands (core/run.R, scripts/ejob), which print summaries.

Pattern-based, not a security boundary: see the cascadia-private skill.
"""
import json
import os
import re
import sys


def deny(reason):
    print(json.dumps({"hookSpecificOutput": {
        "hookEventName": "PreToolUse",
        "permissionDecision": "deny",
        "permissionDecisionReason": reason}}))
    sys.exit(0)


def main():
    try:
        event = json.load(sys.stdin)
    except Exception:
        sys.exit(0)
    proj = os.environ.get("CLAUDE_PROJECT_DIR") or event.get("cwd") or os.getcwd()
    cfg_path = os.path.join(proj, ".cascadia", "private.json")
    if not os.path.isfile(cfg_path):
        sys.exit(0)
    try:
        cfg = json.load(open(cfg_path))
    except Exception:
        cfg = {"private": True, "data_dirs": ["data_raw", "data"]}
    if not cfg.get("private", True):
        sys.exit(0)
    dirs = [d.strip("/").strip() for d in cfg.get("data_dirs", ["data_raw", "data"]) if d.strip()]
    abs_dirs = [os.path.realpath(os.path.join(proj, d)) for d in dirs]
    tool = event.get("tool_name", "")
    ti = event.get("tool_input", {}) or {}
    msg = ("Data-private mode is on for this project ({}): raw records stay out of the conversation. "
           "Use CascadiaTools commands (core/run.R fit / check) that report summaries instead.")

    def inside(path):
        if not path:
            return False
        p = os.path.realpath(path if os.path.isabs(path) else os.path.join(proj, path))
        return any(p == a or p.startswith(a + os.sep) for a in abs_dirs)

    if tool in ("Read", "Edit", "Write", "NotebookEdit"):
        if inside(ti.get("file_path") or ti.get("notebook_path")):
            deny(msg.format(", ".join(dirs)))
    elif tool in ("Grep", "Glob"):
        if inside(ti.get("path") or "") or any(
                re.search(r"(^|/)" + re.escape(d) + r"(/|$)", ti.get("pattern", "") or "") for d in dirs):
            deny(msg.format(", ".join(dirs)))
    elif tool == "Bash":
        cmd = ti.get("command", "")
        mentions = any(re.search(r"(^|[\s'\"=/(])" + re.escape(d) + r"(/|\b)", cmd) for d in dirs) or \
            any(a in cmd for a in abs_dirs)
        pipeline = re.search(r"core/run\.R\s+(fit|figures|simulate|recover|power|selftest|describe)\b", cmd) or \
            re.search(r"scripts/ejob\b.*core/run\.R", cmd)
        if mentions and not pipeline:
            deny(msg.format(", ".join(dirs)))
    sys.exit(0)


if __name__ == "__main__":
    main()
