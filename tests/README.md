# Tests

This directory is intentionally separate from the main analysis workflow.

Current split:

- `qc_result/`
  - stored result snapshots for comparison
- `server_validation/`
  - server-side validation scripts such as Snakemake dry-run checks
- `unit/`
  - lightweight compatibility and boundary checks

Rules:

- Do not wire files under `tests/` into the main workflow configs.
- Do not use test-only absolute paths as defaults in initialization scripts or README examples.
- Run real Snakemake validation on the server side with explicit config paths.
- Keep result snapshots under `tests/` as test artifacts, not as workflow inputs.
