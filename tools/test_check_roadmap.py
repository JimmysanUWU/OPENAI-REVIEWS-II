#!/usr/bin/env python3
"""Negative controls for check_roadmap.py: each corruption of roadmap.json must be reported.

Run from the repository root: python3 tools/test_check_roadmap.py
"""
import copy
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_roadmap as C  # noqa: E402


def errors(data):
    rep = C.Report()
    C.validate(data, rep)
    return rep.errors


def node(data, nid):
    return next(n for n in data["nodes"] if n["id"] == nid)


def main():
    with open(C.JSON_PATH, encoding="utf-8") as fh:
        base = json.load(fh)
    assert errors(base) == [], errors(base)

    corruptions = {
        "unknown dependency": lambda d: node(d, "ETH.I.MAIN")["depends_on"].append("ETH.I.NOPE"),
        "dependency cycle": lambda d: node(d, "ETH.I.SWAP")["depends_on"].append("ETH.I.MAIN"),
        "duplicate id": lambda d: d["nodes"].append(copy.deepcopy(node(d, "KAK.CAP"))),
        "proved on partial": lambda d: node(d, "KAK.T12")["depends_on"].append("KAK.APP.B"),
        "proved on heuristic": lambda d: (node(d, "KAK.CAP").__setitem__("status", ["HEURISTIC"]),
                                          node(d, "KAK.T12")["depends_on"].append("KAK.CAP")),
        "missing file": lambda d: node(d, "KAK.PROP")["where"].append({"file": "reviews/claude/kakeya/99_NOPE.md", "label": "§1"}),
        "wrong theorem label": lambda d: node(d, "ETH.I.DEFCOUNT")["where"].__setitem__(0, {"file": node(d, "ETH.I.DEFCOUNT")["where"][0]["file"], "label": "Thm 4.9"}),
        "wrong section label": lambda d: node(d, "KAK.FIB")["where"].__setitem__(0, {"file": node(d, "KAK.FIB")["where"][0]["file"], "label": "§17"}),
        "missing script": lambda d: node(d, "KAK.CAP")["checks"].append("reviews/claude/kakeya/scripts/nope.py"),
        "FORMAL without name": lambda d: node(d, "KAK.PROP")["status"].append("FORMAL"),
        "unknown Lean name": lambda d: node(d, "KAK.KD").__setitem__("formal", "OAI.Analysis.Fourier.LocalSmoothing.Comparator.kakeya_set"),
        "unknown status": lambda d: node(d, "KAK.CAP")["status"].append("MAYBE"),
        "target without owner": lambda d: d["targets"][0].update(status="IN_PROGRESS", owner=None),
        "done without deliverable": lambda d: d["targets"][0].update(status="DONE"),
        "bad review": lambda d: node(d, "KAK.CAP").__setitem__("reviews", [{"by": "gpt", "verdict": "maybe", "file": "README.md"}]),
        "missing sketch": lambda d: node(d, "KAK.PROP").__setitem__("sketch", "short"),
        "path out of order": lambda d: d["reading_paths"][0]["nodes"].reverse(),
        "node on no path": lambda d: [p["nodes"].remove("KAK.OB.WEAK") for p in d["reading_paths"] if "KAK.OB.WEAK" in p["nodes"]],
        "unknown assumed node": lambda d: d["reading_paths"][3]["assumes"].append("KAK.NOPE"),
        "bad alternative proof": lambda d: node(d, "KAK.T12")["alt_proofs"][0]["depends_on"].append("KAK.NOPE"),
        "alt proof makes a cycle": lambda d: node(d, "KAK.ARROW").__setitem__("alt_proofs", [{"depends_on": ["KAK.T11"], "note": "x"}]),
        "glossary bad pointer": lambda d: d["glossary"][0].__setitem__("label", "§99"),
    }
    failed = []
    for name, corrupt in corruptions.items():
        d = copy.deepcopy(base)
        corrupt(d)
        errs = errors(d)
        status = "detected" if errs else "MISSED"
        print("%-26s %s%s" % (name, status, (": " + errs[0]) if errs else ""))
        if not errs:
            failed.append(name)

    # a stale generated block must be reported, and --write must repair it
    with open(C.MD_PATH, encoding="utf-8") as fh:
        md = fh.read()
    by_id = {n["id"]: n for n in base["nodes"]}
    rendered = C.blocks(base, by_id)
    stale = md.replace("| `KAK.CAP` |", "| `KAK.CAP` | stale", 1)
    rep = C.Report()
    repaired = C.sync(stale, rendered, rep, write=False)
    print("%-26s %s" % ("stale block", "detected" if rep.errors else "MISSED"))
    if not rep.errors:
        failed.append("stale block")
    rep = C.Report()
    if C.sync(stale, rendered, rep, write=True) != md or rep.errors:
        failed.append("repair")
    print("%-26s %s" % ("repair by --write", "ok" if "repair" not in failed else "FAILED"))
    del repaired

    if failed:
        print("FAILED:", failed)
        return 1
    print("ALL OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
