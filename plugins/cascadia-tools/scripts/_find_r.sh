# Sourced by the Rscript and R wrappers: sets R_BIN to the program named in $1.
# Order: CASCADIATOOLS_R (a program or its folder), the lab cluster's R, the PATH, then the
# places installers use. The Windows R installer does not add R to the PATH, and apps opened from
# a Mac's Dock or Finder can have a shorter PATH than Terminal.
find_r() {
  local name="$1" d c
  R_BIN=""
  if [[ -n "$CASCADIATOOLS_R" ]]; then
    if [[ -d "$CASCADIATOOLS_R" ]]; then
      for c in "$CASCADIATOOLS_R/$name" "$CASCADIATOOLS_R/$name.exe" "$CASCADIATOOLS_R/bin/$name" "$CASCADIATOOLS_R/bin/$name.exe"; do
        [[ -f "$c" ]] && { R_BIN="$c"; return; }
      done
      return
    fi
    R_BIN="$CASCADIATOOLS_R"; return
  fi
  [[ -x "/local/cqls/R/R-4.5.3/$name" ]] && { R_BIN="/local/cqls/R/R-4.5.3/$name"; return; }
  R_BIN="$(command -v "$name" 2>/dev/null)"; [[ -n "$R_BIN" ]] && return
  # Windows, in Git Bash
  local roots=() v
  for v in "$PROGRAMFILES" "$ProgramW6432" "${LOCALAPPDATA:+$LOCALAPPDATA/Programs}" "/c/Program Files"; do
    [[ -z "$v" ]] && continue
    command -v cygpath >/dev/null 2>&1 && v="$(cygpath -u "$v")"
    roots+=("$v")
  done
  local found="" nl=$'\n'
  for d in "${roots[@]}"; do
    for c in "$d"/R/R-*/bin/"$name".exe; do
      [[ -f "$c" ]] && found+="$c$nl"
    done
  done
  if [[ -n "$found" ]]; then
    R_BIN="$(printf '%s' "$found" | sort -V | tail -n 1)"; return
  fi
  # macOS
  for c in "/Library/Frameworks/R.framework/Resources/bin/$name" "/usr/local/bin/$name" "/opt/homebrew/bin/$name"; do
    [[ -x "$c" ]] && { R_BIN="$c"; return; }
  done
}
