#!/usr/bin/env python3
"""Validate roadmap.json and keep the generated parts of ROADMAP.md and PROOF_MAP.md in sync.

Usage (from the repository root):
    python3 tools/check_roadmap.py            # validate; exit 1 on any error
    python3 tools/check_roadmap.py --write    # validate, then regenerate the generated blocks
    python3 tools/check_roadmap.py --node ID  # print a node, its proof sketch, and everything it rests on
    python3 tools/check_roadmap.py --path ID  # print a reading path in order (statements and sketches)

Checks:
- schema: required fields, id formats, known kinds and status codes;
- references: every dependency, target base and script exists; every cited file exists;
- labels: structured labels ("§3", "Lemma 2.1", "Thms 4.6", "K3", ...) resolve to a heading or a bold
  statement in the cited file (free-text labels are allowed and not checked);
- logic: the dependency graph (with alternative proofs) is acyclic; a node proved in prose depends only on
  nodes proved in prose; every reading path lists nodes in dependency order (each node's dependencies, for
  its main proof or one alternative, appear earlier or are declared as assumed) and every node is on a path;
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
MAP_PATH = os.path.join(ROOT, "PROOF_MAP.md")

NODE_CODES = {"FORMAL", "PROVED", "PROVED_CITED", "CHECKED", "PARTIAL", "HEURISTIC", "NEW"}
PROSE = {"PROVED", "PROVED_CITED"}
TARGET_CODES = {"OPEN", "IN_PROGRESS", "DONE"}
KINDS = {"theorem", "proposition", "lemma", "corollary", "framework"}
TARGET_KINDS = {"distillation", "research", "formal", "verification"}
OWNERS = {None, "claude", "gpt"}
REVIEW_VERDICTS = {"agree", "dispute", "partial"}
NODE_REQUIRED = ["id", "project", "group", "kind", "title", "statement", "sketch", "status", "where", "depends_on", "checks",
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

def check_where(owner, wheres, rep, texts):
    for w in wheres:
        path = w.get("file", "")
        if not os.path.isfile(os.path.join(ROOT, path)):
            rep.err("%s: missing file %s" % (owner, path))
            continue
        text = texts.setdefault(path, read(path))
        for piece in label_pieces(w.get("label", "")):
            ok = resolve_piece(piece, text)
            if ok is None:
                rep.unchecked += 1
            elif not ok:
                rep.err("%s: label %r not found in %s" % (owner, piece, path))


def proofs(n):
    """Dependency sets of every proof of a node: the main one, then the alternatives."""
    return [n.get("depends_on", [])] + [a.get("depends_on", []) for a in n.get("alt_proofs", [])]


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
        check_where(nid, n.get("where", []), rep, texts)
        if len(n.get("sketch", "")) < 40:
            rep.err("%s: proof sketch missing or too short" % nid)
        for a in n.get("alt_proofs", []):
            if not a.get("depends_on") or not a.get("note"):
                rep.err("%s: alternative proof needs depends_on and note" % nid)
            for d in a.get("depends_on", []):
                if d not in by_id or d == nid:
                    rep.err("%s: bad alternative dependency %s" % (nid, d))
            check_where(nid + " (alt)", a.get("where", []), rep, texts)
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
            for deps in proofs(n):
                for d in deps:
                    if d in by_id and not set(by_id[d]["status"]) & (PROSE | {"FORMAL"}):
                        rep.err("%s is proved in prose but depends on %s (%s)" % (nid, d, ",".join(by_id[d]["status"])))

    # acyclicity (iterative DFS)
    color = {i: 0 for i in by_id}
    for start in by_id:
        if color[start]:
            continue
        stack = [(start, iter(sum(proofs(by_id[start]), [])))]
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
                    stack.append((nxt, iter(sum(proofs(by_id[nxt]), []))))

    # reading paths: dependency order and coverage
    paths = data.get("reading_paths", [])
    pids = [p.get("id") for p in paths]
    if len(set(pids)) != len(pids):
        rep.err("duplicate reading path ids")
    covered = set()
    for p in paths:
        pid = p.get("id", "?")
        if not re.fullmatch(r"P\.(ETH|KAK)\.[A-Z0-9]+", pid):
            rep.err("%s: bad reading path id" % pid)
        known = set(p.get("assumes", []))
        for a in known:
            if a not in by_id:
                rep.err("%s: unknown assumed node %s" % (pid, a))
        for x in p.get("nodes", []):
            if x not in by_id:
                rep.err("%s: unknown node %s" % (pid, x))
                continue
            if not any(set(deps) <= known for deps in proofs(by_id[x])):
                missing = sorted(set(by_id[x]["depends_on"]) - known)
                rep.err("%s: %s comes before its dependencies %s" % (pid, x, missing))
            known.add(x)
            covered.add(x)
    if paths:
        for nid in by_id:
            if nid not in covered:
                rep.err("%s is on no reading path" % nid)
    for g in data.get("glossary", []):
        if not g.get("term") or not g.get("meaning") or g.get("project") not in ("ETH", "KAK"):
            rep.err("bad glossary entry %r" % g.get("term"))
        check_where("glossary %r" % g.get("term"), [{"file": g.get("file", ""), "label": g.get("label", "")}], rep, texts)

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


def prose_escape(s):
    """Escaping for running text (outside tables): keep the source readable as plain text."""
    s = s.replace("*", "\\*")
    return re.sub(r"<(?=[A-Za-z/!])", "&lt;", s)


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
        for a in n.get("alt_proofs", []):
            for d in a["depends_on"]:
                out.append("  %s -.->|alt| %s" % (d.replace(".", "_"), n["id"].replace(".", "_")))
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


BLOCK_RE = re.compile(r"(<!-- BEGIN GENERATED: ([a-z0-9-]+) -->\n)(.*?)(<!-- END GENERATED: \2 -->)", re.S)


def sync(md, rendered, rep, write, fname="ROADMAP.md"):
    found = set()

    def sub(m):
        name = m.group(2)
        found.add(name)
        if name not in rendered:
            rep.err("%s: unknown generated block %r" % (fname, name))
            return m.group(0)
        body = rendered[name] + "\n"
        if m.group(3) != body and not write:
            rep.err("%s: block %r is out of date (run with --write)" % (fname, name))
        return m.group(1) + body + m.group(4)

    new = BLOCK_RE.sub(sub, md)
    for name in rendered:
        if name not in found:
            rep.err("%s: missing generated block %r" % (fname, name))
    return new


def print_node(n):
    print("%s [%s] %s" % (n["id"], ",".join(n["status"]), n["title"]))
    print("    statement: " + n["statement"])
    print("    sketch:    " + n["sketch"])
    print("    uses:      " + (", ".join(n["depends_on"]) or "-"))
    for a in n.get("alt_proofs", []):
        print("    alt proof: uses " + ", ".join(a["depends_on"]) + " -- " + a["note"])
    print("    where:     " + "; ".join("%s (%s)" % (w["file"], w["label"]) for w in n["where"]))


def show(nid, by_id):
    if nid not in by_id:
        print("unknown node", nid)
        return 1
    for x in [nid] + sorted(ancestors(nid, by_id)):
        print_node(by_id[x])
    print("prose chain:", chain(nid, by_id))
    return 0


def show_path(pid, data, by_id):
    p = next((q for q in data.get("reading_paths", []) if q["id"] == pid), None)
    if p is None:
        print("unknown path", pid, "- known:", ", ".join(q["id"] for q in data.get("reading_paths", [])))
        return 1
    print("%s: %s\n%s\nassumes: %s\n" % (p["id"], p["title"], p["purpose"], ", ".join(p["assumes"]) or "-"))
    for i, x in enumerate(p["nodes"], 1):
        print("%d. " % i, end="")
        print_node(by_id[x])
    return 0


def anchor(pid, nid):
    return (pid + "--" + nid).lower().replace(".", "-")


def canonical_anchors(data):
    out = {}
    for p in data.get("reading_paths", []):
        for x in p["nodes"]:
            out.setdefault(x, anchor(p["id"], x))
    return out


def link_md(path):
    return "[%s](%s)" % (short(path), path)


def render_path(p, by_id, canon):
    here = set(p["nodes"])
    known = set(p["assumes"])
    lines = ["*%s*" % prose_escape(p["purpose"])]
    if p["assumes"]:
        lines.append("")
        lines.append("Assumes: " + ", ".join("[`%s`](#%s)" % (a, canon[a]) for a in p["assumes"]) + ".")
    for i, x in enumerate(p["nodes"], 1):
        n = by_id[x]
        def ref(d):
            return "[`%s`](#%s)" % (d, anchor(p["id"], d) if d in here else canon[d])
        lines += ["", '#### <a id="%s"></a>%d. `%s` — %s' % (anchor(p["id"], x), i, x, prose_escape(n["title"])), ""]
        meta = "*%s*" % ", ".join(n["status"])
        ch = chain(x, by_id)
        if ch != "—":
            meta += " · prose chain: %s" % ch
        if n["depends_on"]:
            meta += " · uses " + ", ".join(ref(d) for d in n["depends_on"])
        lines.append(meta)
        lines += ["", "**Statement.** " + prose_escape(n["statement"]), "", "**Proof sketch.** " + prose_escape(n["sketch"])]
        if not set(n["depends_on"]) <= known:
            lines += ["", "**On this path the result is reached by the alternative proof below**; the sketch above is the main proof."]
        known.add(x)
        for a in n.get("alt_proofs", []):
            lines += ["", "**Alternative proof.** Uses " + ", ".join(ref(d) for d in a["depends_on"]) + ". " + prose_escape(a["note"])
                      + (" (" + "; ".join("%s %s" % (link_md(w["file"]), prose_escape(w["label"])) for w in a.get("where", [])) + ")" if a.get("where") else "")]
        tail = "**Full proof:** " + "; ".join("%s %s" % (link_md(w["file"]), prose_escape(w["label"])) for w in n["where"])
        if n["checks"]:
            tail += " · **checks:** " + ", ".join(link_md(c) for c in n["checks"])
        if n["formal"]:
            tail += " · **Lean:** `%s`" % n["formal"]
        if n.get("notes"):
            tail += " · *Note:* " + prose_escape(n["notes"])
        lines += ["", tail]
    return "\n".join(lines)


def render_path_index(data):
    rows = ["| Path | Read it for | Steps | Assumes |", "|---|---|---|---|"]
    for p in data.get("reading_paths", []):
        rows.append("| [%s](#%s) | %s | %d | %s |" % (
            md_escape(p["title"]), p["id"].lower().replace(".", "-"), prose_escape(p["purpose"]).replace("|", "\\|"), len(p["nodes"]),
            ", ".join("`%s`" % a for a in p["assumes"]) or "—"))
    return "\n".join(rows)


def render_glossary(data, project):
    rows = ["| Term | Meaning | Defined in |", "|---|---|---|"]
    for g in data.get("glossary", []):
        if g["project"] == project:
            rows.append("| %s | %s | %s %s |" % (md_escape(g["term"]), md_escape(g["meaning"]), link_md(g["file"]), md_escape(g["label"])))
    return "\n".join(rows)


def render_node_index(data, by_id):
    canon = canonical_anchors(data)
    onpaths = {}
    for p in data.get("reading_paths", []):
        for x in p["nodes"]:
            onpaths.setdefault(x, []).append(p["id"])
    rows = ["| Node | Title | Status | On paths |", "|---|---|---|---|"]
    for n in data["nodes"]:
        rows.append("| [`%s`](#%s) | %s | %s | %s |" % (n["id"], canon[n["id"]], md_escape(n["title"]), ", ".join(n["status"]),
                                                        ", ".join(onpaths.get(n["id"], []))))
    return "\n".join(rows)


def map_blocks(data, by_id):
    canon = canonical_anchors(data)
    out = {"path-index": render_path_index(data), "glossary-kak": render_glossary(data, "KAK"),
           "glossary-eth": render_glossary(data, "ETH"), "node-index": render_node_index(data, by_id)}
    for p in data.get("reading_paths", []):
        out["path-" + p["id"].lower().replace(".", "-")] = render_path(p, by_id, canon)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--write", action="store_true", help="regenerate the generated blocks of ROADMAP.md")
    ap.add_argument("--node", help="print a node and its transitive dependencies")
    ap.add_argument("--path", help="print a reading path in order")
    args = ap.parse_args()
    with open(JSON_PATH, encoding="utf-8") as fh:
        data = json.load(fh)
    rep = Report()
    by_id = validate(data, rep)
    if args.node:
        return show(args.node, by_id)
    if args.path:
        return show_path(args.path, data, by_id)
    if rep.errors:
        print("\n".join("ERROR " + e for e in rep.errors))
        print("FAILED: %d error(s); ROADMAP.md not touched" % len(rep.errors))
        return 1
    for path, rendered in ((MD_PATH, blocks(data, by_id)), (MAP_PATH, map_blocks(data, by_id))):
        with open(path, encoding="utf-8") as fh:
            md = fh.read()
        new = sync(md, rendered, rep, args.write, os.path.basename(path))
        if args.write and new != md:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(new)
            print("%s regenerated" % os.path.basename(path))
    if rep.errors:
        print("\n".join("ERROR " + e for e in rep.errors))
        print("FAILED: %d error(s)" % len(rep.errors))
        return 1
    print("OK: %d nodes, %d targets, %d free-text labels not checked"
          % (len(data["nodes"]), len(data["targets"]), rep.unchecked))
    return 0


if __name__ == "__main__":
    sys.exit(main())
