#!/usr/bin/env python3
"""Validate roadmap.json and keep the generated parts of ROADMAP.md in sync.

Usage (from the repository root):
    python3 tools/check_roadmap.py            # validate; exit 1 on any error
    python3 tools/check_roadmap.py --write    # validate, then regenerate the generated blocks of ROADMAP.md
    python3 tools/check_roadmap.py --node ID  # print a node and everything it rests on

Checks:
- schema: required fields, id formats, known kinds and status codes;
- references: every dependency, target base and script exists; every cited file exists;
- labels: structured labels ("§3", "Lemma 2.1", "Thms 4.6", "K3", ...) resolve to a heading or a bold
  statement in the cited file (free-text labels are allowed and not checked);
- logic: the dependency graph is acyclic; a node proved in prose depends only on nodes proved in prose;
  FORMAL marks exactly the nodes with an audited Lean name, and that name occurs in the pinned statement file;
- sync: the blocks of ROADMAP.md between <!-- BEGIN GENERATED: x --> and <!-- END GENERATED: x -->
  equal what this script renders from roadmap.json.
Standard library only; deterministic.
"""
import argparse
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JSON_PATH = os.path.join(ROOT, "roadmap.json")
MD_PATH = os.path.join(ROOT, "ROADMAP.md")

NODE_CODES = {"FORMAL", "PROVED", "PROVED_CITED", "CHECKED", "PARTIAL", "HEURISTIC", "NEW"}
PROSE = {"PROVED", "PROVED_CITED"}
TARGET_CODES = {"OPEN", "IN_PROGRESS", "DONE"}
KINDS = {"theorem", "proposition", "lemma", "corollary", "framework"}
TARGET_KINDS = {"distillation", "research", "formal", "verification"}
OWNERS = {None, "claude", "gpt"}
REVIEW_VERDICTS = {"agree", "dispute", "partial"}
NODE_REQUIRED = ["id", "project", "group", "kind", "title", "statement", "status", "where", "depends_on", "checks",
                 "paper_ref", "formal", "confidence"]
TARGET_REQUIRED = ["id", "kind", "title", "goal", "builds_on", "status", "owner", "difficulty", "acceptance"]
FORMAL_SOURCES = {"ETH": "source/eth/SubexponentialSAT.statement.lean",
                  "KAK": "source/kakeya/LocalSmoothing.challenge.lean"}
WORD = {"Lemma": "Lemma", "Lemmas": "Lemma", "Thm": "Theorem", "Thms": "Theorem", "Theorem": "Theorem",
        "Theorems": "Theorem", "Prop.": "Proposition", "Props.": "Proposition", "Proposition": "Proposition",
        "Cor.": "Corollary", "Cors.": "Corollary", "Corollary": "Corollary", "Remark": "Remark",
        "Definition": "Definition"}


class Report:
    def __init__(self):
        self.errors = []
        self.unchecked = 0

    def err(self, msg):
        self.errors.append(msg)


def read(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as fh:
        return fh.read()


# ---------------------------------------------------------------- labels

def label_pieces(label):
    """Split 'X, Y (comment), Z' into pieces, dropping trailing parentheticals."""
    out = []
    for piece in label.split(","):
        piece = re.sub(r"\s*\([^()]*\)\s*$", "", piece.strip())
        if piece:
            out.append(piece)
    return out


def resolve_piece(piece, text):
    """Return True/False if the piece is structured and found/missing, None if free text."""
    m = re.fullmatch(r"§(\d+)(?:\.(\d+|[A-Z]))?(?:-§?(\d+))?", piece)
    if m:
        secs = [m.group(1)] if m.group(3) is None else [m.group(1), m.group(3)]
        if m.group(2) is not None:
            pat = r"^#{2,4} %s\.%s[ .]" % (m.group(1), m.group(2))
            return re.search(pat, text, re.M) is not None
        return all(re.search(r"^#{2,4} %s\. " % s, text, re.M) for s in secs)
    m = re.fullmatch(r"(Lemmas?|Thms?|Theorems?|Props?\.|Proposition|Cors?\.|Corollary|Remark|Definition)\s+"
                     r"(\d+\.\d+)(?:-(\d+\.\d+))?", piece)
    if m:
        word = WORD[m.group(1)]
        nums = [m.group(2)] + ([m.group(3)] if m.group(3) else [])
        return all(re.search(r"\*\*%s %s\b" % (word, re.escape(x)), text) for x in nums)
    m = re.fullmatch(r"([EK]\d)", piece)
    if m:
        return re.search(r"^## %s\. " % m.group(1), text, re.M) is not None
    m = re.fullmatch(r"Thm (I|II|III|IV|V)", piece)
    if m:
        return re.search(r"Theorem %s\b" % m.group(1), text) is not None
    return None


# ---------------------------------------------------------------- validation

def validate(data, rep):
    for key in ["schema", "updated", "status_codes", "projects", "nodes", "targets"]:
        if key not in data:
            rep.err("missing top-level key %r" % key)
    nodes, targets = data.get("nodes", []), data.get("targets", [])
    ids = [n.get("id") for n in nodes]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        rep.err("duplicate node ids: %s" % sorted(dup))
    by_id = {n["id"]: n for n in nodes if "id" in n}
    texts = {}
    formal_texts = {p: read(f) for p, f in FORMAL_SOURCES.items()}

    for n in nodes:
        nid = n.get("id", "?")
        for f in NODE_REQUIRED:
            if f not in n:
                rep.err("%s: missing field %r" % (nid, f))
        if not re.fullmatch(r"(ETH|KAK)(\.[A-Z0-9]+)+", nid):
            rep.err("%s: bad id" % nid)
        if n.get("project") != nid.split(".")[0]:
            rep.err("%s: project does not match id" % nid)
        if n.get("kind") not in KINDS:
            rep.err("%s: unknown kind %r" % (nid, n.get("kind")))
        st = n.get("status", [])
        if not st or not set(st) <= NODE_CODES:
            rep.err("%s: bad status %r" % (nid, st))
        for d in n.get("depends_on", []):
            if d not in by_id:
                rep.err("%s: unknown dependency %s" % (nid, d))
            if d == nid:
                rep.err("%s: depends on itself" % nid)
        for w in n.get("where", []):
            path = w.get("file", "")
            if not os.path.isfile(os.path.join(ROOT, path)):
                rep.err("%s: missing file %s" % (nid, path))
                continue
            text = texts.setdefault(path, read(path))
            for piece in label_pieces(w.get("label", "")):
                ok = resolve_piece(piece, text)
                if ok is None:
                    rep.unchecked += 1
                elif not ok:
                    rep.err("%s: label %r not found in %s" % (nid, piece, path))
        for c in n.get("checks", []):
            if not os.path.isfile(os.path.join(ROOT, c)):
                rep.err("%s: missing check script %s" % (nid, c))
        formal = n.get("formal")
        if ("FORMAL" in st) != (formal is not None):
            rep.err("%s: FORMAL status and formal name must go together" % nid)
        if formal:
            ns, _, name = formal.rpartition(".")
            src = formal_texts[n["project"]]
            if not re.search(r"^theorem %s\b" % re.escape(name), src, re.M) or ("namespace " + ns) not in src:
                rep.err("%s: formal name %s not found in %s" % (nid, formal, FORMAL_SOURCES[n["project"]]))
        for r in n.get("reviews", []):
            if r.get("by") not in {"claude", "gpt"} or r.get("verdict") not in REVIEW_VERDICTS:
                rep.err("%s: bad review entry %r" % (nid, r))
            elif not os.path.isfile(os.path.join(ROOT, r.get("file", ""))):
                rep.err("%s: review file missing %s" % (nid, r.get("file")))
        conf = n.get("confidence")
        if conf is not None and not (0 <= conf <= 1):
            rep.err("%s: confidence outside [0,1]" % nid)
        if set(st) & PROSE:
            for d in n.get("depends_on", []):
                if d in by_id and not set(by_id[d]["status"]) & (PROSE | {"FORMAL"}):
                    rep.err("%s is proved in prose but depends on %s (%s)" % (nid, d, ",".join(by_id[d]["status"])))

    # acyclicity (iterative DFS)
    color = {i: 0 for i in by_id}
    for start in by_id:
        if color[start]:
            continue
        stack = [(start, iter(by_id[start].get("depends_on", [])))]
        color[start] = 1
        while stack:
            node, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                color[node] = 2
                stack.pop()
            elif nxt in color:
                if color[nxt] == 1:
                    rep.err("dependency cycle through %s -> %s" % (node, nxt))
                elif color[nxt] == 0:
                    color[nxt] = 1
                    stack.append((nxt, iter(by_id[nxt].get("depends_on", []))))

    tids = [t.get("id") for t in targets]
    if len(set(tids)) != len(tids):
        rep.err("duplicate target ids")
    for t in targets:
        tid = t.get("id", "?")
        for f in TARGET_REQUIRED:
            if f not in t:
                rep.err("%s: missing field %r" % (tid, f))
        if not re.fullmatch(r"T\.(ETH|KAK|X)\.\d+", tid):
            rep.err("%s: bad target id" % tid)
        if t.get("kind") not in TARGET_KINDS:
            rep.err("%s: unknown kind %r" % (tid, t.get("kind")))
        if t.get("status") not in TARGET_CODES:
            rep.err("%s: bad status %r" % (tid, t.get("status")))
        if t.get("owner") not in OWNERS:
            rep.err("%s: unknown owner %r" % (tid, t.get("owner")))
        if t.get("status") == "IN_PROGRESS" and not t.get("owner"):
            rep.err("%s: in progress without an owner" % tid)
        for b in t.get("builds_on", []):
            if b not in by_id:
                rep.err("%s: unknown base node %s" % (tid, b))
        if t.get("status") == "DONE":
            dl = t.get("deliverable")
            if not dl or not os.path.isfile(os.path.join(ROOT, dl)):
                rep.err("%s: DONE but deliverable missing" % tid)
    return by_id


# ---------------------------------------------------------------- rendering

def ancestors(nid, by_id):
    seen, todo = set(), [nid]
    while todo:
        x = todo.pop()
        for d in by_id[x].get("depends_on", []):
            if d not in seen:
                seen.add(d)
                todo.append(d)
    return seen


def chain(nid, by_id):
    """Weakest prose link among the node and its ancestors."""
    st = set(by_id[nid]["status"])
    if not st & PROSE:
        return "—"
    cited = sorted(a for a in ancestors(nid, by_id) | {nid} if "PROVED_CITED" in by_id[a]["status"])
    return "closed" if not cited else "cites " + ", ".join(cited)


def md_escape(s):
    return s.replace("|", "\\|").replace("*", "\\*").replace("<", "&lt;").replace(">", "&gt;")


def short(path):
    m = re.match(r"reviews/claude/(eth|kakeya)/(?:scripts/)?([^/]+)$", path)
    if not m:
        return path
    name = m.group(2)
    if name.endswith(".py"):
        return name[:-3]
    return "%s/%s" % (m.group(1), name.split("_")[0])


def render_nodes(project, by_id, nodes):
    rows = ["| ID | Result | Status | Prose chain | Depends on | Where | Checks | Lean |",
            "|---|---|---|---|---|---|---|---|"]
    for n in nodes:
        if n["project"] != project:
            continue
        where = "; ".join("[%s](%s) %s" % (short(w["file"]), w["file"], md_escape(w["label"])) for w in n["where"])
        checks = ", ".join("[%s](%s)" % (short(c), c) for c in n["checks"]) or "—"
        lean = ("`%s`" % n["formal"].rpartition(".")[2]) if n["formal"] else "—"
        result = "**%s.** %s" % (md_escape(n["title"]), md_escape(n["statement"]))
        if n.get("notes"):
            result += " *Note:* " + md_escape(n["notes"])
        status = ", ".join(n["status"])
        for r in n.get("reviews", []):
            status += "; [%s: %s](%s)" % (r["by"], r["verdict"], r["file"])
        rows.append("| `%s` | %s | %s | %s | %s | %s | %s | %s |" % (
            n["id"], result, status, chain(n["id"], by_id),
            ", ".join("`%s`" % d for d in n["depends_on"]) or "—", where, checks, lean))
    return "\n".join(rows)


def mermaid_label(s):
    s = s.replace('"', "'").replace("{", "").replace("}", "")
    return "".join({"#": "#35;", "<": "#60;", ">": "#62;"}.get(ch, ch) for ch in s)


def mermaid_class(n):
    st = set(n["status"])
    if "FORMAL" in st:
        return "formal"
    if "HEURISTIC" in st:
        return "heuristic"
    if "PARTIAL" in st:
        return "partial"
    if "PROVED_CITED" in st:
        return "cited"
    return "proved"


def render_graph(project, nodes):
    mine = [n for n in nodes if n["project"] == project]
    groups = []
    for n in mine:
        if n["group"] not in groups:
            groups.append(n["group"])
    out = ["```mermaid", "flowchart TD"]
    for gi, g in enumerate(groups):
        out.append('  subgraph %s_g%d["%s"]' % (project, gi, mermaid_label(g)))
        for n in mine:
            if n["group"] == g:
                new = " (new)" if "NEW" in n["status"] else ""
                out.append('    %s["%s<br/>%s%s"]' % (n["id"].replace(".", "_"), n["id"], mermaid_label(n["title"]), new))
        out.append("  end")
    for n in mine:
        for d in n["depends_on"]:
            out.append("  %s --> %s" % (d.replace(".", "_"), n["id"].replace(".", "_")))
    for cls in ["formal", "proved", "cited", "partial", "heuristic"]:
        members = [n["id"].replace(".", "_") for n in mine if mermaid_class(n) == cls]
        if members:
            out.append("  class %s %s" % (",".join(members), cls))
    out += ["  classDef formal fill:#c8e6c9,stroke:#2e7d32,color:#000",
            "  classDef proved fill:#bbdefb,stroke:#1565c0,color:#000",
            "  classDef cited fill:#fff9c4,stroke:#f9a825,color:#000",
            "  classDef partial fill:#ffccbc,stroke:#d84315,color:#000",
            "  classDef heuristic fill:#e1bee7,stroke:#6a1b9a,color:#000,stroke-dasharray:4 3",
            "```"]
    return "\n".join(out)


def render_targets(targets):
    rows = ["| ID | Kind | Title | Status | Owner | Builds on | Goal | Acceptance |", "|---|---|---|---|---|---|---|---|"]
    for t in targets:
        goal = md_escape(t["goal"])
        if t.get("deliverable"):
            goal += " Deliverable: [`%s`](%s)." % (t["deliverable"].rpartition("/")[2], t["deliverable"])
        rows.append("| `%s` | %s | %s | %s | %s | %s | %s | %s |" % (
            t["id"], t["kind"], md_escape(t["title"]), t["status"], t["owner"] or "—",
            ", ".join("`%s`" % b for b in t["builds_on"]), goal, "; ".join(md_escape(a) for a in t["acceptance"])))
    return "\n".join(rows)


def render_summary(nodes, targets):
    rows = ["| Project | Results | FORMAL | PROVED | PROVED_CITED | PARTIAL | HEURISTIC | NEW | CHECKED |",
            "|---|---|---|---|---|---|---|---|---|"]
    for p in ["ETH", "KAK"]:
        mine = [n for n in nodes if n["project"] == p]
        cnt = lambda c: sum(c in n["status"] for n in mine)
        rows.append("| %s | %d | %d | %d | %d | %d | %d | %d | %d |" % (
            p, len(mine), cnt("FORMAL"), cnt("PROVED"), cnt("PROVED_CITED"), cnt("PARTIAL"), cnt("HEURISTIC"),
            cnt("NEW"), cnt("CHECKED")))
    rows.append("")
    rows.append("Targets: %d open, %d in progress, %d done." % tuple(
        sum(t["status"] == s for t in targets) for s in ["OPEN", "IN_PROGRESS", "DONE"]))
    return "\n".join(rows)


def blocks(data, by_id):
    nodes, targets = data["nodes"], data["targets"]
    return {"summary": render_summary(nodes, targets),
            "eth-graph": render_graph("ETH", nodes), "eth-nodes": render_nodes("ETH", by_id, nodes),
            "kak-graph": render_graph("KAK", nodes), "kak-nodes": render_nodes("KAK", by_id, nodes),
            "targets": render_targets(targets)}


BLOCK_RE = re.compile(r"(<!-- BEGIN GENERATED: ([a-z-]+) -->\n)(.*?)(<!-- END GENERATED: \2 -->)", re.S)


def sync(md, rendered, rep, write):
    found = set()

    def sub(m):
        name = m.group(2)
        found.add(name)
        if name not in rendered:
            rep.err("ROADMAP.md: unknown generated block %r" % name)
            return m.group(0)
        body = rendered[name] + "\n"
        if m.group(3) != body and not write:
            rep.err("ROADMAP.md: block %r is out of date (run with --write)" % name)
        return m.group(1) + body + m.group(4)

    new = BLOCK_RE.sub(sub, md)
    for name in rendered:
        if name not in found:
            rep.err("ROADMAP.md: missing generated block %r" % name)
    return new


def show(nid, by_id):
    if nid not in by_id:
        print("unknown node", nid)
        return 1
    for x in [nid] + sorted(ancestors(nid, by_id)):
        n = by_id[x]
        print("%s [%s] %s" % (x, ",".join(n["status"]), n["title"]))
        print("    " + n["statement"])
        print("    where: " + "; ".join("%s (%s)" % (w["file"], w["label"]) for w in n["where"]))
    print("prose chain:", chain(nid, by_id))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--write", action="store_true", help="regenerate the generated blocks of ROADMAP.md")
    ap.add_argument("--node", help="print a node and its transitive dependencies")
    args = ap.parse_args()
    with open(JSON_PATH, encoding="utf-8") as fh:
        data = json.load(fh)
    rep = Report()
    by_id = validate(data, rep)
    if args.node:
        return show(args.node, by_id)
    if rep.errors:
        print("\n".join("ERROR " + e for e in rep.errors))
        print("FAILED: %d error(s); ROADMAP.md not touched" % len(rep.errors))
        return 1
    with open(MD_PATH, encoding="utf-8") as fh:
        md = fh.read()
    new = sync(md, blocks(data, by_id), rep, args.write)
    if args.write and new != md:
        with open(MD_PATH, "w", encoding="utf-8") as fh:
            fh.write(new)
        print("ROADMAP.md regenerated")
    if rep.errors:
        print("\n".join("ERROR " + e for e in rep.errors))
        print("FAILED: %d error(s)" % len(rep.errors))
        return 1
    print("OK: %d nodes, %d targets, %d free-text labels not checked"
          % (len(data["nodes"]), len(data["targets"]), rep.unchecked))
    return 0


if __name__ == "__main__":
    sys.exit(main())
