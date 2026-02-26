#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

PATCH_ROOT = Path("/home/user/Codex/patches")


def _now_tag() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _slugify(s: str, maxlen: int = 48) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    s = re.sub(r"-+", "-", s).strip("-")
    return s[:maxlen] or "patch"


def _run(cmd: List[str], cwd: Optional[Path] = None) -> str:
    p = subprocess.run(cmd, cwd=str(cwd) if cwd else None, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nstdout:\n{p.stdout}\nstderr:\n{p.stderr}")
    return p.stdout


def repo_root() -> Path:
    out = _run(["git", "rev-parse", "--show-toplevel"]).strip()
    return Path(out)


def current_commit() -> str:
    return _run(["git", "rev-parse", "HEAD"]).strip()


def diff_patch(paths: Optional[List[str]] = None) -> str:
    cmd = ["git", "diff"]
    if paths:
        cmd += ["--"] + paths
    return _run(cmd)


def ensure_clean_index() -> None:
    # disallow staged changes for safety; working tree may be dirty (we’ll diff it)
    out = _run(["git", "diff", "--cached", "--name-only"]).strip()
    if out:
        raise RuntimeError("Staged changes detected. Please commit or unstage before creating patch.")


@dataclass
class PatchMeta:
    skill: str
    reason: str
    base_commit: str
    timestamp: str
    files: List[str]
    run_ref: Optional[str] = None
    evidence: Optional[Dict] = None


def create_patch(
    skill: str,
    reason: str,
    files: Optional[List[str]] = None,
    run_ref: Optional[str] = None,
    evidence: Optional[Dict] = None,
) -> Path:
    ensure_clean_index()

    root = repo_root()
    base = current_commit()
    ts = _now_tag()
    slug = _slugify(reason)

    patch_dir = PATCH_ROOT / skill / f"{ts}_{slug}"
    patch_dir.mkdir(parents=True, exist_ok=False)

    # compute patch
    patch_text = diff_patch(files)
    if not patch_text.strip():
        # no changes => no patch
        shutil.rmtree(patch_dir)
        raise RuntimeError("No diff found. Nothing to propose as a patch.")

    (patch_dir / "diff.patch").write_text(patch_text, encoding="utf-8")

    meta = PatchMeta(
        skill=skill,
        reason=reason,
        base_commit=base,
        timestamp=datetime.now().isoformat(),
        files=files or [],
        run_ref=run_ref,
        evidence=evidence,
    )
    (patch_dir / "meta.json").write_text(json.dumps(meta.__dict__, indent=2), encoding="utf-8")

    # changelog
    lines = []
    lines.append(f"# Patch Proposal: {skill}")
    lines.append("")
    lines.append(f"- timestamp: {meta.timestamp}")
    lines.append(f"- base_commit: {meta.base_commit}")
    lines.append(f"- reason: {meta.reason}")
    if meta.run_ref:
        lines.append(f"- run_ref: {meta.run_ref}")
    if files:
        lines.append("- files_changed:")
        for f in files:
            lines.append(f"  - {f}")
    else:
        lines.append("- files_changed: (git diff determines actual changes)")
    if evidence:
        lines.append("- evidence:")
        # pretty json block
        ev = json.dumps(evidence, indent=2, ensure_ascii=False)
        lines.append("```json")
        lines.append(ev)
        lines.append("```")
    lines.append("")
    lines.append("## How to apply")
    lines.append("```bash")
    lines.append(f"python tools/patch_manager.py apply --dir {patch_dir}")
    lines.append("```")
    (patch_dir / "CHANGELOG.md").write_text("\n".join(lines), encoding="utf-8")

    return patch_dir


def apply_patch(patch_dir: Path) -> None:
    patch_file = patch_dir / "diff.patch"
    if not patch_file.exists():
        raise FileNotFoundError(f"diff.patch not found in {patch_dir}")

    # Safety: ensure we are on the same base commit (optional but recommended)
    meta_path = patch_dir / "meta.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        base = meta.get("base_commit")
        if base and base != current_commit():
            raise RuntimeError(
                f"Base commit mismatch.\n"
                f"Patch base: {base}\nCurrent: {current_commit()}\n"
                f"Checkout the patch base commit or cherry-pick manually."
            )

    # Apply
    _run(["git", "apply", "--index", str(patch_file)])
    print(f"[OK] Applied patch to index: {patch_file}")
    print("Next: commit the changes:")
    print('  git commit -m "Apply patch proposal: ..."')
    print("  git push")


def main():
    ap = argparse.ArgumentParser(prog="patch_manager.py")
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("create")
    c.add_argument("--skill", required=True)
    c.add_argument("--reason", required=True)
    c.add_argument("--files", default=None, help="comma-separated paths to include in diff (optional)")
    c.add_argument("--run_ref", default=None)
    c.add_argument("--evidence_json", default=None, help="path to evidence json file (optional)")
    c.set_defaults(which="create")

    a = sub.add_parser("apply")
    a.add_argument("--dir", required=True, help="patch directory containing diff.patch/meta.json")
    a.set_defaults(which="apply")

    args = ap.parse_args()

    if args.which == "create":
        files = [x.strip() for x in args.files.split(",")] if args.files else None
        evidence = None
        if args.evidence_json:
            evidence = json.loads(Path(args.evidence_json).read_text())
        p = create_patch(
            skill=args.skill,
            reason=args.reason,
            files=files,
            run_ref=args.run_ref,
            evidence=evidence,
        )
        print(f"[OK] Created patch: {p}")
        print(f"- changelog: {p/'CHANGELOG.md'}")
        print(f"- diff: {p/'diff.patch'}")
        print(f"- meta: {p/'meta.json'}")

    elif args.which == "apply":
        apply_patch(Path(args.dir))


if __name__ == "__main__":
    main()
