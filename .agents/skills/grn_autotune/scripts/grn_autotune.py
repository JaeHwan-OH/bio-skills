#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, os, time
from pathlib import Path
import yaml

def load_paramspace(p: Path) -> dict:
    return yaml.safe_load(p.read_text()) or {}

def load_tf_expr(export_dir: Path):
    # tf_expr_by_group.tsv: gene, group, avg_expr, pct_expr
    rows = []
    with (export_dir / "tf_expr_by_group.tsv").open() as f:
        header = f.readline().strip().split("\t")
        idx = {k:i for i,k in enumerate(header)}
        for line in f:
            x = line.rstrip("\n").split("\t")
            rows.append((x[idx["gene"]], x[idx["group"]], float(x[idx["avg_expr"]]), float(x[idx["pct_expr"]])))
    return rows

def product_dict(d: dict):
    # simple cartesian over two lists (gate only)
    gate = d.get("tf_expression_gate", {})
    pct_list = gate.get("min_pct_expr", [0.05])
    avg_list = gate.get("min_avg_expr", [0.10])
    for pct in pct_list:
        for avg in avg_list:
            yield {"min_pct_expr": float(pct), "min_avg_expr": float(avg)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--export_dir", required=True)
    ap.add_argument("--paramspace", default=str(Path(__file__).resolve().parents[1] / "templates" / "paramspace.default.yaml"))
    ap.add_argument("--out_dir", default=None)
    args = ap.parse_args()

    export_dir = Path(args.export_dir)
    if args.out_dir:
        out_dir = Path(args.out_dir)
    else:
        out_dir = export_dir.parent.parent / ("autotune_" + time.strftime("%Y%m%d_%H%M%S"))
    results = out_dir / "results"
    repro = out_dir / "repro"
    results.mkdir(parents=True, exist_ok=True)
    repro.mkdir(parents=True, exist_ok=True)

    ps = load_paramspace(Path(args.paramspace))
    rows = load_tf_expr(export_dir)

    leaderboard = []
    for cfg in product_dict(ps):
        min_pct = cfg["min_pct_expr"]
        min_avg = cfg["min_avg_expr"]
        # count TFs passing gate in ANY group (first minimal metric)
        passed = set()
        for gene, group, avg, pct in rows:
            if pct >= min_pct and avg >= min_avg:
                passed.add(gene)
        score = len(passed)
        leaderboard.append({"min_pct_expr": min_pct, "min_avg_expr": min_avg, "score_tf_pass": score})

    leaderboard.sort(key=lambda x: x["score_tf_pass"], reverse=True)
    best = leaderboard[0] if leaderboard else {}

    # write leaderboard
    with (results / "leaderboard.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["min_pct_expr","min_avg_expr","score_tf_pass"])
        w.writeheader()
        for r in leaderboard:
            w.writerow(r)

    (results / "best_params.yaml").write_text(yaml.safe_dump({"tf_expression_gate": best}, sort_keys=False))
    (repro / "tuning_report.md").write_text(
        "# grn-autotune report\n\n"
        f"- export_dir: {export_dir}\n"
        f"- paramspace: {args.paramspace}\n\n"
        "현재 버전은 TF 발현 게이트 튜닝의 최소 뼈대입니다.\n"
        "다음 단계에서 footprint/binding + GRN 스코어(다운스트림 품질 지표)를 추가해, 실제 성능 최적화를 수행합니다.\n"
    )
    print("[OK]", out_dir)

if __name__ == "__main__":
    main()
