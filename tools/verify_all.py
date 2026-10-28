#!/usr/bin/env python3
"""Run every check script in the repository and report one line per script.

Usage (from the repository root):
    python3 tools/verify_all.py              # everything (about 3 minutes with 4 jobs)
    python3 tools/verify_all.py --quick      # skip the scripts that take more than a minute
    python3 tools/verify_all.py --project KAK
    python3 tools/verify_all.py --jobs 1     # sequential

A script passes if it exits with status 0 and prints "ALL OK" or "ALL PASS". The two report-only scripts
(verify_cost_floor, verify_parameters) compute tables rather than pass/fail checks; they pass on exit
status 0. The roadmap checker and its negative controls are run first. Exit status 1 if anything fails.
Requirements: Python 3.10+; the Kakeya scripts also need numpy, sympy, scipy, mpmath.
"""
import argparse
import concurrent.futures as cf
import glob
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLOW = {"verify_cost_floor.py", "verify_counting_pipeline.py", "verify_wave_current.py"}
REPORT_ONLY = {"verify_cost_floor.py", "verify_parameters.py"}
DIRS = {"ETH": "reviews/claude/eth/scripts", "KAK": "reviews/claude/kakeya/scripts"}


def nodes_checked_by():
    with open(os.path.join(ROOT, "roadmap.json"), encoding="utf-8") as fh:
        data = json.load(fh)
    out = {}
    for n in data["nodes"]:
        for c in n["checks"]:
            out.setdefault(os.path.basename(c), []).append(n["id"])
    return out


def run(cmd, cwd, timeout):
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        out = p.stdout + p.stderr
        code = p.returncode
    except subprocess.TimeoutExpired:
        out, code = "TIMEOUT", -1
    return code, out, time.time() - t0


def job(path, timeout):
    name = os.path.basename(path)
    code, out, dt = run([sys.executable, name], os.path.dirname(path), timeout)
    ok = code == 0 and (name in REPORT_ONLY or "ALL OK" in out or "ALL PASS" in out)
    last = [l for l in out.strip().splitlines() if l.strip()][-1:] or [""]
    return name, ok, dt, last[0][:100]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="skip scripts slower than a minute")
    ap.add_argument("--project", choices=["ETH", "KAK"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=900)
    args = ap.parse_args()

    failures = 0
    for cmd in (["tools/check_roadmap.py"], ["tools/test_check_roadmap.py"]):
        code, out, dt = run([sys.executable] + cmd, ROOT, 120)
        ok = code == 0
        failures += not ok
        print("%-4s %-34s %6.1fs  %s" % ("ok" if ok else "FAIL", cmd[0], dt, out.strip().splitlines()[-1][:100]))

    scripts = []
    for proj, d in DIRS.items():
        if args.project and proj != args.project:
            continue
        for path in sorted(glob.glob(os.path.join(ROOT, d, "verify_*.py"))):
            if args.quick and os.path.basename(path) in SLOW:
                continue
            scripts.append((proj, path))
    checked = nodes_checked_by()
    with cf.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(job, path, args.timeout): proj for proj, path in scripts}
        results = []
        for f in cf.as_completed(futs):
            results.append((futs[f],) + f.result())
    for proj, name, ok, dt, last in sorted(results):
        failures += not ok
        nodes = ", ".join(checked.get(name, [])) or "(toolkit/companion)"
        print("%-4s %-3s %-34s %6.1fs  %s" % ("ok" if ok else "FAIL", proj, name, dt, last))
        print("         checks: %s" % nodes)
    print("ALL OK" if not failures else "FAILED: %d" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
