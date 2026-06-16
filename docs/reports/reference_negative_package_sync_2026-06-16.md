# Reference Negative Package Sync 2026-06-16

## Scope

This report records the local synchronization status for the negative reference FASTQ package requested for downstream reference-modeling work.

## Source

- Remote host: `fengxian`
- Remote package: `/data/project/CNV/PGT-A/rawdata/lib_test/reference/reference_negative_samples_pass_warn_20260615.tar`
- Remote MD5 file: `/data/project/CNV/PGT-A/rawdata/lib_test/reference/reference_negative_samples_pass_warn_20260615.tar.md5`
- Package summary:
  - sample count: `32`
  - FASTQ file count: `64`
  - XX count: `20`
  - XY count: `12`
  - FASTQ total: `96.886 GiB`

## Local Target

- Local directory: `D:\Pipeline\PGT-A\reference_packages`
- Local tar: `D:\Pipeline\PGT-A\reference_packages\reference_negative_samples_pass_warn_20260615.tar`
- Local remote-MD5 copy: `D:\Pipeline\PGT-A\reference_packages\reference_negative_samples_pass_warn_20260615.tar.md5`
- Local MD5 result: `D:\Pipeline\PGT-A\reference_packages\reference_negative_samples_pass_warn_20260615.local_md5.txt`
- Local verification JSON: `D:\Pipeline\PGT-A\reference_packages\reference_negative_samples_pass_warn_20260615.md5_verify.json`

## Verification

Completed on local storage after the remote transfer finished.

| check | value |
| --- | --- |
| expected bytes | `104030259200` |
| actual bytes | `104030259200` |
| expected MD5 | `fa7036739ecec4685364734762d55d08` |
| actual MD5 | `fa7036739ecec4685364734762d55d08` |
| state | `verified` |
| checked at | `2026-06-16T06:42:11` |

## Notes

- The local package is intentionally excluded from Git via `reference_packages/` in `.gitignore`.
- This report is the committed audit record; the tar, PID files, transfer logs, and local MD5 sidecars remain local data artifacts.
