#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

RUN_ROOT = Path("/home/user/Codex/runs")
REGISTRY = Path("/home/user/Codex/_registry/registry.json")


def now_tag() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def load_registry() -> dict:
    ensure_dir(REGISTRY.parent)
    if REGISTRY.exists():
        return json.loads(REGISTRY.read_text())
    return {"datasets": {}, "fragments_seen": [], "rmd_seen": []}


def save_registry(reg: dict) -> None:
    ensure_dir(REGISTRY.parent)
    REGISTRY.write_text(json.dumps(reg, indent=2, sort_keys=True))


def read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="ignore")


# -------------------------
# Rmd mining (paths)
# -------------------------
def extract_paths_from_rmd(rmd_path: Path) -> Dict[str, List[str]]:
    """
    Extract likely absolute paths from Rmd.
    """
    text = read_text(rmd_path)

    # quoted absolute paths
    candidates = re.findall(r'["\'](/[^"\']+)["\']', text)

    frags = [c for c in candidates if re.search(r"fragments\.tsv\.gz(\.bgz)?$", c)]
    beds = [c for c in candidates if c.endswith(".bed")]
    h5 = [c for c in candidates if c.endswith((".h5", ".h5ad", ".h5mu"))]
    outs = [c for c in candidates if ("/outs" in c or "cellranger" in c.lower())]

    # unquoted fragments patterns
    frags += re.findall(r'(/[^ \n\t"]+fragments\.tsv\.gz(?:\.bgz)?)', text)

    def uniq(xs):
        seen = set()
        out = []
        for x in xs:
            x = x.strip()
            if not x or x in seen:
                continue
            seen.add(x)
            out.append(x)
        return out

    return {
        "fragments": uniq(frags),
        "beds": uniq(beds),
        "h5": uniq(h5),
        "outs": uniq(outs),
        "all_candidates": uniq(candidates),
    }


def pick_existing(paths: List[str], must_have_tbi: bool = False) -> Optional[str]:
    for p in paths:
        pp = Path(p)
        if pp.exists():
            if must_have_tbi:
                if Path(str(pp) + ".tbi").exists():
                    return str(pp)
            else:
                return str(pp)
    return None


# -------------------------
# h5ad introspection
# -------------------------
@dataclass
class ColGuess:
    stage_col: str
    genotype_col: str
    cluster_col: str
    stage_value: str
    genotype_value: str


def _values(adata, col: str, limit: int = 2000) -> List[str]:
    vals = adata.obs[col].astype(str).values
    if len(vals) > limit:
        vals = vals[:limit]
    return list(vals)


def _score_stage(values: List[str], target: str) -> Tuple[int, Optional[str]]:
    t = target.lower()
    uniq = sorted(set(v.lower() for v in values))
    # exact
    if t in uniq:
        return 100, target
    # substring
    for u in uniq:
        if t in u:
            return 80, u
    # embryonic pattern
    for u in uniq:
        if re.match(r"^e\d+(\.\d+)?$", u):
            # if same day number matches
            if u.replace("e", "").startswith(target.lower().replace("e", "")):
                return 60, u
    return 0, None


def _score_genotype(values: List[str], target: str) -> Tuple[int, Optional[str]]:
    t = target.lower()
    uniq = sorted(set(v.lower() for v in values))
    if t in uniq:
        return 100, target
    # common synonyms
    syn = []
    if t == "wt":
        syn = ["wildtype", "wild-type", "control", "ctrl", "w/t"]
    if t in ("ko", "cko"):
        syn = ["knockout", "cKO", "mut", "mutant"]
    for s in syn:
        for u in uniq:
            if s.lower() == u or s.lower() in u:
                return 70, u
    # partial match (WT_KO labels)
    for u in uniq:
        if t in u:
            return 60, u
    return 0, None


def guess_columns(adata, target_stage: str, target_genotype: str) -> ColGuess:
    cols = list(adata.obs.columns)

    # candidate lists
    stage_name_hint = ["stage", "Stage", "embryo_stage", "time", "age", "development", "embryonic"]
    geno_name_hint = ["genotype", "Genotype", "geno", "condition"]
    cluster_name_hint = ["active_ident", "seurat_clusters", "cluster", "clusters", "leiden", "louvain", "celltype"]

    def name_score(c: str, hints: List[str]) -> int:
        cl = c.lower()
        s = 0
        for h in hints:
            if h.lower() == cl:
                s += 50
            elif h.lower() in cl:
                s += 15
        return s

    # Stage
    best_stage = (0, None, None)  # score, col, matched_value
    for c in cols:
        ns = name_score(c, stage_name_hint)
        vs = _values(adata, c)
        cs, mv = _score_stage(vs, target_stage)
        score = ns + cs
        if score > best_stage[0]:
            best_stage = (score, c, mv)

    # Genotype
    best_geno = (0, None, None)
    for c in cols:
        ns = name_score(c, geno_name_hint)
        vs = _values(adata, c)
        cs, mv = _score_genotype(vs, target_genotype)
        score = ns + cs
        if score > best_geno[0]:
            best_geno = (score, c, mv)

    # Cluster (no target value needed)
    best_cluster = (0, None)
    for c in cols:
        ns = name_score(c, cluster_name_hint)
        # heuristics: not too many unique, not too few
        uniq = len(set(_values(adata, c, limit=5000)))
        cs = 0
        if 2 <= uniq <= 200:
            cs = 40
        elif 200 < uniq <= 2000:
            cs = 10
        score = ns + cs
        if score > best_cluster[0]:
            best_cluster = (score, c)

    # Fail if not found
    if not best_stage[1] or not best_stage[2]:
        raise RuntimeError(
            f"Could not auto-detect stage column/value for target '{target_stage}'. "
            f"Please provide --stage_col and/or --stage_value."
        )
    if not best_geno[1] or not best_geno[2]:
        raise RuntimeError(
            f"Could not auto-detect genotype column/value for target '{target_genotype}'. "
            f"Please provide --genotype_col and/or --genotype_value."
        )
    if not best_cluster[1]:
        raise RuntimeError(
            "Could not auto-detect cluster column. Please provide --cluster_col."
        )

    return ColGuess(
        stage_col=best_stage[1],
        genotype_col=best_geno[1],
        cluster_col=best_cluster[1],
        stage_value=str(best_stage[2]),
        genotype_value=str(best_geno[2]),
    )


def write_barcode_lists(adata, cluster_col: str, outdir: Path) -> Dict[str, str]:
    ensure_dir(outdir)
    clusters = sorted(set(adata.obs[cluster_col].astype(str).values))
    mapping = {}
    for cl in clusters:
        mask = adata.obs[cluster_col].astype(str) == cl
        bcs = adata.obs_names[mask].tolist()
        fn = outdir / f"cluster_{cl}.txt"
        fn.write_text("\n".join(bcs) + "\n")
        mapping[str(cl)] = str(fn)
    return mapping


def main():
    ap = argparse.ArgumentParser(prog="multiome_bridge.py")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prep_h5ad", help="Auto-detect columns, subset E10.5_WT (or user-specified), create cluster barcode lists, mine Rmd for paths.")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--rmd", required=True)
    p.add_argument("--project", default="Gdf6_project")
    p.add_argument("--target_stage", default="E10.5")
    p.add_argument("--target_genotype", default="WT")
    p.add_argument("--stage_col", default=None)
    p.add_argument("--genotype_col", default=None)
    p.add_argument("--cluster_col", default=None)
    p.add_argument("--stage_value", default=None)
    p.add_argument("--genotype_value", default=None)
    p.set_defaults(which="prep_h5ad")

    args = ap.parse_args()

    ensure_dir(RUN_ROOT)
    reg = load_registry()

    if args.which == "prep_h5ad":
        h5ad = Path(args.h5ad).expanduser().resolve()
        rmd = Path(args.rmd).expanduser().resolve()
        if not h5ad.exists():
            raise FileNotFoundError(f"h5ad not found: {h5ad}")
        if not rmd.exists():
            raise FileNotFoundError(f"Rmd not found: {rmd}")

        # Run directory
        run_dir = RUN_ROOT / args.project / f"{now_tag()}_prep_h5ad"
        export_dir = run_dir / "02_export"
        logs_dir = run_dir / "logs"
        ensure_dir(export_dir)
        ensure_dir(logs_dir)

        # Load AnnData
        try:
            import anndata as ad
        except Exception as e:
            raise RuntimeError("Missing python package: anndata. Install it in your base/venv.") from e

        adata = ad.read_h5ad(str(h5ad))

        # Auto-detect columns unless user provided
        if args.stage_col and args.genotype_col and args.cluster_col and args.stage_value and args.genotype_value:
            guess = ColGuess(
                stage_col=args.stage_col,
                genotype_col=args.genotype_col,
                cluster_col=args.cluster_col,
                stage_value=args.stage_value,
                genotype_value=args.genotype_value,
            )
        else:
            guess = guess_columns(adata, args.target_stage, args.target_genotype)
            # allow partial overrides
            if args.stage_col: guess.stage_col = args.stage_col
            if args.genotype_col: guess.genotype_col = args.genotype_col
            if args.cluster_col: guess.cluster_col = args.cluster_col
            if args.stage_value: guess.stage_value = args.stage_value
            if args.genotype_value: guess.genotype_value = args.genotype_value

        # Subset
        stage_mask = adata.obs[guess.stage_col].astype(str).str.lower() == str(guess.stage_value).lower()
        geno_mask = adata.obs[guess.genotype_col].astype(str).str.lower().str.contains(str(guess.genotype_value).lower())
        sub = adata[stage_mask & geno_mask].copy()

        if sub.n_obs == 0:
            # write a debug report then fail
            dbg = {
                "h5ad": str(h5ad),
                "stage_col": guess.stage_col,
                "genotype_col": guess.genotype_col,
                "cluster_col": guess.cluster_col,
                "stage_value": guess.stage_value,
                "genotype_value": guess.genotype_value,
                "hint": "subset is empty; check stage/genotype values in adata.obs",
                "obs_columns": list(adata.obs.columns)[:80],
            }
            (run_dir / "DEBUG_empty_subset.json").write_text(json.dumps(dbg, indent=2))
            raise RuntimeError("Subset is empty. See DEBUG_empty_subset.json and override --stage_value/--genotype_value or columns.")

        # Save subset h5ad
        subset_path = export_dir / f"{args.target_stage}_{args.target_genotype}.h5ad"
        sub.write_h5ad(str(subset_path))

        # Write barcode lists by cluster
        bc_dir = export_dir / "barcodes_by_cluster"
        bc_map = write_barcode_lists(sub, guess.cluster_col, bc_dir)

        # Mine Rmd for path hints and validate existence
        mined = extract_paths_from_rmd(rmd)
        picked_frag = pick_existing(mined["fragments"], must_have_tbi=True)
        picked_bed = pick_existing(mined["beds"], must_have_tbi=False)

        validation = {
            "fragments_candidates": mined["fragments"],
            "fragments_picked": picked_frag,
            "peaks_bed_candidates": mined["beds"],
            "peaks_bed_picked": picked_bed,
            "h5_candidates": mined["h5"],
            "outs_candidates": mined["outs"],
        }
        (export_dir / "rmd_mined_paths.json").write_text(json.dumps(validation, indent=2))

        # registry update
        ds_key = f"{h5ad.stem}__{args.target_stage}_{args.target_genotype}"
        reg["datasets"][ds_key] = {
            "h5ad": str(h5ad),
            "rmd": str(rmd),
            "subset_h5ad": str(subset_path),
            "stage_col": guess.stage_col,
            "genotype_col": guess.genotype_col,
            "cluster_col": guess.cluster_col,
            "stage_value": guess.stage_value,
            "genotype_value": guess.genotype_value,
            "fragments_hint": picked_frag,
            "peaks_bed_hint": picked_bed,
            "updated": datetime.now().isoformat(),
        }
        if str(rmd) not in reg.get("rmd_seen", []):
            reg.setdefault("rmd_seen", []).append(str(rmd))
        if picked_frag and picked_frag not in reg.get("fragments_seen", []):
            reg.setdefault("fragments_seen", []).append(picked_frag)
        save_registry(reg)

        # Write a human-readable summary
        summary = {
            "run_dir": str(run_dir),
            "subset_h5ad": str(subset_path),
            "n_obs": int(sub.n_obs),
            "n_vars": int(sub.n_vars),
            "detected": guess.__dict__,
            "barcodes_by_cluster_dir": str(bc_dir),
            "num_clusters": len(bc_map),
            "rmd_paths": validation,
            "dataset_key": ds_key,
        }
        (run_dir / "SUMMARY.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
