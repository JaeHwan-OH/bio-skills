---
name: multiome-bridge
description: >
  Seurat/Signac multiome 통합 rds(이미 QC/클러스터링 완료)에서 구조를 자동 탐색해
  scPrinter 또는 SCENIC+ 실행에 필요한 번들을 export하고, conda 없이 venv(+pip)로
  tool별 환경을 캐시(/media/user/HDD1/data)하여 CPU-only로 실행까지 자동화한다.
---

# 고정 경로 규칙
- RUN_ROOT: /home/user/Codex/runs
- REGISTRY: /home/user/Codex/_registry/registry.json
- CACHE_ROOT: /media/user/HDD1/data
- VENV_ROOT: /media/user/HDD1/data/_envs
- PIP_CACHE: /media/user/HDD1/data/_cache/pip
- SCPRINTER_DATA: /media/user/HDD1/data/_cache/scprinter

# 불변 원칙
1) 추측 금지: rds/fragments 경로가 없으면 요구 후 중단
2) Preflight 우선: 파일 존재/권한/python3.10/3.11 존재 확인
3) 환경 재사용: venv는 (tool + python버전 + requirements hash)로 캐시 재사용
4) 재현성 저장: repro/commands.sh, repro/versions.txt, repro/run_export.R, repro/run_tool.py
5) 통합 object 대응: fragments가 여러 개일 수 있음 (barcode -> fragments map)
6) group 라벨: active.ident(Idents) 사용
7) TF 발현 게이트용 통계 export 포함:
   - tf_expr_by_group.tsv (avg_expr, pct_expr)
