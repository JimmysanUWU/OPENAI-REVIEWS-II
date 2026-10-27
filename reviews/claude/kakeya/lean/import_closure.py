"""Compute the OAI-internal import closure of a Lean module in a blobless clone.

Fetches missing blobs level by level with `git checkout HEAD -- <paths>`, then scans
every file in the closure for proof escape hatches.
Usage: python3 -I import_closure.py <repo_root> <root.module> [<root.module> ...]
"""
import json
import os
import re
import subprocess
import sys

repo = sys.argv[1]
roots = sys.argv[2:]
lean = os.path.join(repo, "lean")
tree = set(subprocess.run(["git", "-C", repo, "ls-tree", "-r", "--name-only", "HEAD", "lean/OAI"],
                          capture_output=True, text=True, check=True).stdout.split())


def path_of(mod):
    return "lean/" + mod.replace(".", "/") + ".lean"


IMPORT = re.compile(r"^\s*import\s+(.+)$")
seen, frontier, external = set(), list(roots), set()
while frontier:
    need = [path_of(m) for m in frontier if not os.path.exists(os.path.join(repo, path_of(m)))]
    need = [p for p in need if p in tree]
    for i in range(0, len(need), 400):
        subprocess.run(["git", "-C", repo, "checkout", "HEAD", "--"] + need[i:i + 400],
                       check=True, capture_output=True)
    nxt = []
    for m in frontier:
        if m in seen:
            continue
        seen.add(m)
        p = os.path.join(repo, path_of(m))
        if not os.path.exists(p):
            external.add(m)
            continue
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                mm = IMPORT.match(line)
                if mm:
                    for dep in mm.group(1).split():
                        if dep.startswith("OAI.") and dep not in seen:
                            nxt.append(dep)
                        elif not dep.startswith("OAI."):
                            external.add(dep.split(".")[0])
    frontier = nxt

files = sorted(path_of(m) for m in seen if os.path.exists(os.path.join(repo, path_of(m))))
patterns = {
    "sorry": r"\bsorry\b",
    "admit": r"\badmit\b",
    "axiom_decl": r"^\s*(private\s+|protected\s+)?axiom\s",
    "native_decide": r"\bnative_decide\b",
    "implemented_by": r"implemented_by",
    "extern": r"@\[extern",
    "unsafe": r"\bunsafe\b",
    "ofReduceBool": r"ofReduceBool",
    "skipKernelTC": r"skipKernelTC",
    "opaque_decl": r"^\s*opaque\s",
    "set_option_debug": r"set_option\s+debug\.",
}
hits = {k: [] for k in patterns}
lines = 0
for f in files:
    with open(os.path.join(repo, f), encoding="utf-8") as fh:
        text = fh.read()
    lines += text.count("\n")
    # strip comments for the scan
    code = re.sub(r"/-.*?-/", "", text, flags=re.S)
    code = re.sub(r"--[^\n]*", "", code)
    for k, pat in patterns.items():
        for mm in re.finditer(pat, code, flags=re.M):
            ln = code[:mm.start()].count("\n") + 1
            hits[k].append(f"{f}:{ln}")
print(json.dumps({"roots": roots, "oai_modules": len(files), "lines": lines,
                  "external_roots": sorted(external - set(seen)),
                  "hits": {k: v[:20] for k, v in hits.items()},
                  "hit_counts": {k: len(v) for k, v in hits.items()}}, indent=1))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "closure_files.txt"), "w") as fh:
    fh.write("\n".join(files))
