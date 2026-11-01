"""Run the original checks with a forward-slash __file__ on Windows.

The mathematics and original scripts are unchanged. Full outputs are retained.
"""
from pathlib import Path
import concurrent.futures
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'reviews/gpt/check-logs'
OUT.mkdir(parents=True, exist_ok=True)
WRAPPER = "import runpy,sys; p=sys.argv[1]; sys.path.insert(0,p.rsplit('/',1)[0]); sys.argv=[p]; runpy.run_path(p,run_name='__main__')"

def run(path):
    t0 = time.monotonic()
    proc = subprocess.run([sys.executable, '-c', WRAPPER, path.as_posix()],
                          cwd=path.parent, capture_output=True, text=True,
                          encoding='utf-8', timeout=900)
    output = proc.stdout + proc.stderr
    (OUT / (path.stem + '.txt')).write_text(output, encoding='utf-8')
    report = path.name in {'verify_parameters.py', 'verify_cost_floor.py'}
    lines = [s.strip() for s in proc.stdout.splitlines() if s.strip()]
    ok = proc.returncode == 0 and (report or lines[-1:] in [['ALL OK'], ['ALL PASS']])
    result = dict(script=path.relative_to(ROOT).as_posix(), exit_code=proc.returncode,
                  ok=ok, seconds=round(time.monotonic()-t0,2), last=lines[-1:] or [''])
    print(json.dumps(result), flush=True)
    return result

if __name__ == '__main__':
    paths = sorted((ROOT/'reviews/claude/eth/scripts').glob('verify_*.py'))
    paths += sorted((ROOT/'reviews/claude/kakeya/scripts').glob('verify_*.py'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(run, paths))
    (OUT/'summary.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    sys.exit(0 if all(r['ok'] for r in results) else 1)
