"""Validate the explicitly documented API without importing scientific backends."""
from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
API_DIR = Path(__file__).resolve().parent
SECTIONS = re.compile(r"(?m)^([A-Za-z][A-Za-z ]+)\n-{3,}\s*$")


def sections(doc: str) -> dict[str, str]:
    matches = list(SECTIONS.finditer(doc))
    return {
        match.group(1): doc[match.end():matches[i + 1].start() if i + 1 < len(matches) else len(doc)]
        for i, match in enumerate(matches)
    }


def locate(qualified: str) -> ast.FunctionDef | ast.ClassDef:
    parts = qualified.split(".")
    if len(parts) < 3 or parts[0] != "pygrnwang":
        raise ValueError("unsupported qualified name")
    module = ROOT / "pygrnwang" / (parts[1] + ".py")
    node = ast.parse(module.read_text(encoding="utf-8-sig"), filename=str(module))
    for name in parts[2:]:
        matches = [
            item for item in node.body
            if isinstance(item, (ast.FunctionDef, ast.ClassDef)) and item.name == name
        ]
        if len(matches) != 1:
            raise ValueError(f"cannot uniquely find {name}")
        node = matches[0]
    return node


def parameters(node: ast.FunctionDef | ast.ClassDef) -> set[str]:
    if isinstance(node, ast.ClassDef):
        constructors = [
            item for item in node.body
            if isinstance(item, ast.FunctionDef) and item.name == "__init__"
        ]
        if not constructors:
            return set()
        node = constructors[0]
    args = node.args
    names = {a.arg for a in args.posonlyargs + args.args + args.kwonlyargs}
    names.update(a.arg for a in (args.vararg, args.kwarg) if a is not None)
    return names - {"self", "cls"}


def main() -> int:
    manifest = json.loads((API_DIR / "public-api.json").read_text(encoding="utf-8"))
    errors = []
    seen = set()
    for group, entries in manifest["groups"].items():
        page = (API_DIR / (group + ".md")).read_text(encoding="utf-8")
        for qualified in entries:
            if qualified in seen:
                errors.append(f"{qualified}: duplicate manifest entry")
            seen.add(qualified)
            try:
                node = locate(qualified)
            except (OSError, SyntaxError, ValueError) as exc:
                errors.append(f"{qualified}: {exc}")
                continue
            doc = ast.get_docstring(node) or ""
            parsed = sections(doc)
            text = parsed.get("Parameters", "") + "\n" + parsed.get("Other Parameters", "")
            documented = {
                name.strip().lstrip("*")
                for line in text.splitlines()
                if re.match(r"^[*\w][*\w, ]*\s*:\s*\S", line)
                for name in line.split(":", 1)[0].split(",")
            }
            missing = parameters(node) - documented
            if missing:
                errors.append(f"{qualified}: undocumented arguments: {', '.join(sorted(missing))}")
            if not any(parsed.get(section, "").strip() for section in ("Returns", "Yields")):
                errors.append(f"{qualified}: missing Returns or Yields")
            if not doc.split("\n", 1)[0].strip():
                errors.append(f"{qualified}: missing summary")
            if qualified not in page:
                errors.append(f"{qualified}: absent from {group}.md")
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"API coverage passed: {len(seen)} explicit objects; all signature parameters and returns documented.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
