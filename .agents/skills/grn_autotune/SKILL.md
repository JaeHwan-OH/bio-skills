---
name: grn-autotune
description: >
  export bundle(tf_expr_by_group 포함)를 입력으로 받아 TF 발현 게이트를 강제 적용하고,
  downstream(footprint/GRN) 파라미터 탐색(sweep)의 실행/집계/리포팅을 자동화한다.
---

# 기본 게이트 (변경 가능)
- min_pct_expr: 0.05
- min_avg_expr: 0.10  (Seurat data slot 평균)

# 입력 요구
export_dir 안에 최소:
- tf_expr_by_group.tsv (gene, group, avg_expr, pct_expr)
- cells.tsv (barcode, group)

# 산출물
- results/leaderboard.csv
- results/best_params.yaml
- repro/tuning_report.md
