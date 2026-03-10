#!/usr/bin/env python3
import argparse
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def run_text(command, cwd=None):
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        output = (result.stdout or "").strip()
        if not output:
            output = f"<empty> (exit={result.returncode})"
        return output.replace("\n", " | ")
    except Exception as exc:
        return f"<error: {exc}>"


def main():
    parser = argparse.ArgumentParser(description="Collect reproducibility metadata for pipeline runs.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--fastp", required=True)
    parser.add_argument("--bwa", required=True)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--wisecondorx", required=True)
    parser.add_argument("--python-bin", required=True)
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    project_root = Path(args.project_root)

    rows = []
    rows.append(("generated_utc", datetime.now(timezone.utc).isoformat()))
    rows.append(("project_root", str(project_root.resolve())))
    rows.append(("git_branch", run_text(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=project_root)))
    rows.append(("git_commit", run_text(["git", "rev-parse", "HEAD"], cwd=project_root)))
    rows.append(("git_status_short", run_text(["git", "status", "--short"], cwd=project_root)))
    rows.append(("fastp_version", run_text([args.fastp, "--version"])))
    rows.append(("bwa_version", run_text([args.bwa], cwd=project_root)))
    rows.append(("samtools_version", run_text([args.samtools, "--version"])))
    rows.append(("wisecondorx_version", run_text([args.wisecondorx, "--version"])))
    rows.append(("python_version", run_text([args.python_bin, "--version"])))

    with open(output, "w", encoding="utf-8") as handle:
        handle.write("key\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    main()
