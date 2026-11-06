# CascadiaTools data-private guard: called by private_guard, with the hook's JSON event on stdin.
# Denies tool calls that would show raw records from the project's listed data folders: file tools
# on those paths, and shell commands that mention them, except CascadiaTools pipeline commands
# (core/run.R, scripts/ejob), which print summaries.
# Pattern-based, not a security boundary: see the cascadia-private skill.
extra <- strsplit(Sys.getenv("CASCADIATOOLS_RLIB", ""), .Platform$path.sep, fixed = TRUE)[[1]]
extra <- c(extra, file.path(Sys.getenv("CASCADIATOOLS_ROOT", ""), "Rlib"))
.libPaths(c(extra[nzchar(extra) & dir.exists(extra)], .libPaths()))
if (!requireNamespace("jsonlite", quietly = TRUE)) quit(status = 3, save = "no")

is_win <- .Platform$OS.type == "windows"
deny <- function(reason) {
  cat(jsonlite::toJSON(list(hookSpecificOutput = list(
    hookEventName = "PreToolUse", permissionDecision = "deny", permissionDecisionReason = reason)),
    auto_unbox = TRUE))
  quit(status = 0, save = "no")
}
is_abs <- function(p) grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", p)
# Resolve a path that may not exist yet: resolve its deepest existing folder, then re-attach the rest.
real <- function(p) {
  p <- gsub("\\\\", "/", p)
  tail <- character()
  while (!file.exists(p) && dirname(p) != p) {
    tail <- c(basename(p), tail)
    p <- dirname(p)
  }
  out <- do.call(file.path, as.list(c(normalizePath(p, winslash = "/", mustWork = FALSE), tail)))
  if (is_win) tolower(out) else out
}

event <- tryCatch(jsonlite::fromJSON(paste(readLines(file("stdin"), warn = FALSE), collapse = "\n"),
                                     simplifyVector = FALSE),
                  error = function(e) NULL)
if (is.null(event)) quit(status = 3, save = "no")  # refused by private_guard: an unreadable event is not safe to allow
proj <- Sys.getenv("CLAUDE_PROJECT_DIR", "")
if (!nzchar(proj)) proj <- if (is.null(event$cwd)) getwd() else event$cwd
cfg <- tryCatch(jsonlite::fromJSON(file.path(proj, ".cascadia", "private.json"), simplifyVector = FALSE),
                error = function(e) list(private = TRUE, data_dirs = list("data_raw", "data")))
if (identical(cfg$private, FALSE)) quit(status = 0, save = "no")
dirs <- if (is.null(cfg$data_dirs)) c("data_raw", "data") else unlist(cfg$data_dirs)
dirs <- trimws(gsub("^[/\\\\]+|[/\\\\]+$", "", dirs))
dirs <- dirs[nzchar(dirs)]
abs_dirs <- vapply(dirs, function(d) real(file.path(proj, d)), "")
msg <- sprintf(paste("Data-private mode is on for this project (%s): raw records stay out of the conversation.",
                     "Use CascadiaTools commands (core/run.R fit / check) that report summaries instead."),
               paste(dirs, collapse = ", "))
inside <- function(path) {
  if (is.null(path) || !nzchar(path)) return(FALSE)
  p <- real(if (is_abs(path)) path else file.path(proj, path))
  any(p == abs_dirs | startsWith(p, paste0(abs_dirs, "/")))
}
rx <- function(d) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", d)
tool <- if (is.null(event$tool_name)) "" else event$tool_name
ti <- if (is.null(event$tool_input)) list() else event$tool_input

if (tool %in% c("Read", "Edit", "Write", "NotebookEdit")) {
  fp <- if (!is.null(ti$file_path)) ti$file_path else ti$notebook_path
  if (inside(fp)) deny(msg)
} else if (tool %in% c("Grep", "Glob")) {
  pat <- if (is.null(ti$pattern)) "" else ti$pattern
  in_pat <- any(vapply(dirs, function(d) grepl(paste0("(^|[/\\\\])", rx(d), "([/\\\\]|$)"), pat, perl = TRUE), NA))
  if (inside(ti$path) || in_pat) deny(msg)
} else if (tool %in% c("Bash", "PowerShell")) {
  cmd <- if (is.null(ti$command)) "" else ti$command
  mentions <- any(vapply(dirs, function(d)
    grepl(paste0("(^|[\\s'\"=/\\\\(])", rx(d), "([/\\\\]|\\b)"), cmd, perl = TRUE), NA)) ||
    any(vapply(abs_dirs, function(a) grepl(a, gsub("\\\\", "/", cmd), fixed = TRUE, ignore.case = FALSE), NA))
  pipeline <- grepl("core[/\\\\]run\\.R\\s+(fit|figures|simulate|recover|power|selftest|describe)\\b", cmd, perl = TRUE) ||
    grepl("scripts[/\\\\]ejob\\b.*core[/\\\\]run\\.R", cmd, perl = TRUE)
  if (mentions && !pipeline) deny(msg)
}
quit(status = 0, save = "no")
