#!/bin/bash
# Build the remaining ETH modules one at a time (dependencies first) to avoid the memory limit.
# Usage: build_sequential.sh <deadline-seconds>
export PATH=$HOME/.elan/bin:$PATH
S=/tmp/claude-0/-home-user/585f675e-483c-5f5d-b157-c342ce21c632/scratchpad
P=$S/leanbuild/proj
LOG=$S/leanbuild/eth_seq.log
deadline=$(( $(date +%s) + $1 ))
cd "$P" || exit 1
built() { [ -e "$P/.lake/build/lib/lean/$(echo "$1" | tr '.' '/').olean" ]; }
ready() {   # all OAI imports of module $1 are built
  local f="$P/$(echo "$1" | tr '.' '/').lean"
  for dep in $(grep -E '^import OAI\.' "$f" | awk '{print $2}'); do built "$dep" || return 1; done
  return 0
}
while [ "$(date +%s)" -lt "$deadline" ]; do
  progress=0
  for m in $(cat "$S/eth_missing.txt"); do
    built "$m" && continue
    ready "$m" || continue
    left=$(( deadline - $(date +%s) ))
    [ "$left" -lt 30 ] && break 2
    echo "== $m ($(date +%T))" >> "$LOG"
    timeout "$left" lake build "$m" >> "$LOG" 2>&1
    echo "rc=$?" >> "$LOG"
    find .lake/build/ir -name '*.setup.json' -delete
    progress=1
  done
  all=1; for m in $(cat "$S/eth_missing.txt"); do built "$m" || all=0; done
  [ "$all" = 1 ] && { echo "ALL BUILT" >> "$LOG"; break; }
  [ "$progress" = 0 ] && { echo "NO PROGRESS" >> "$LOG"; break; }
done
for m in $(cat "$S/eth_missing.txt"); do built "$m" && echo "built $m" || echo "missing $m"; done
