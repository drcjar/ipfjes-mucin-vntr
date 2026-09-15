"""render_lint.py: flag markdown that pandoc renders differently from how it reads.

Reads pandoc's JSON AST on stdin (pandoc -f markdown -t json draft.md | python3 render_lint.py)
and prints one line per finding; no output means clean. Used by check_paper.sh.

- Emphasis that begins or ends with punctuation/space, or contains a superscript digit:
  the signature of a literal * pairing with a later * (e.g. "Schwartz⁵*, Wain³⁶*" became
  italic ", Wain³⁶" and both markers vanished). Escape literal asterisks as \\*.
- A paragraph joining separate source lines: consecutive lines with no blank line or
  trailing backslash merge (the affiliations ran together). End such lines with \\.
- A citation element: without --citeproc, @text is parsed as a citation. Escape as \\@.
"""
import json
import sys

SUPERSCRIPT = set("⁰¹²³⁴⁵⁶⁷⁸⁹")


def text(x):
    if isinstance(x, dict):
        t, c = x.get("t"), x.get("c")
        if t == "Str":
            return c
        if t in ("Space", "SoftBreak", "LineBreak"):
            return " "
        if t == "Code":
            return c[1]
        return text(c)
    if isinstance(x, list):
        return "".join(text(i) for i in x)
    return ""


findings = []


def walk(x):
    if isinstance(x, dict):
        t, c = x.get("t"), x.get("c")
        if t in ("Emph", "Strong"):
            s = text(c)
            if s[:1] in " ,.;:)" or s[-1:] in " ,;(" or SUPERSCRIPT & set(s):
                findings.append(f"[EMPH] {t} span looks like paired literal asterisks: '{s[:70]}'")
        elif t == "Cite":
            findings.append(f"[CITE] parsed as a citation: '{text(c[1])[:70]}' (escape @ as \\@)")
        elif t in ("Para", "Plain") and any(isinstance(i, dict) and i.get("t") == "SoftBreak" for i in c):
            findings.append(f"[JOINED] separate source lines merged into one paragraph: '{text(c)[:70]}'")
        walk(c)
    elif isinstance(x, list):
        for i in x:
            walk(i)


walk(json.load(sys.stdin)["blocks"])
if findings:
    print("\n".join(findings))
