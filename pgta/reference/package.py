#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Build a frozen reference package.")
    parser.add_argument("--package-name", required=True)
    parser.add_argument("--package-dir", required=True)
    parser.add_argument("--reference-files", nargs="+", required=True)
    parser.add_argument("--asset-files", nargs="+", required=True)
    parser.add_argument("--metadata-files", nargs="*", default=[])
    parser.add_argument("--manifest-output", required=True)
    parser.add_argument("--inventory-output", required=True)
    parser.add_argument("--readme-output", required=True)
    parser.add_argument("--done-output", required=True)
    return parser.parse_args()


def sha256_file(path_value):
    digest = hashlib.sha256()
    with open(path_value, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_group(paths, destination_root, group_name):
    records = []
    group_dir = destination_root / group_name
    group_dir.mkdir(parents=True, exist_ok=True)
    for path_str in paths:
        source = Path(path_str)
        if not source.exists():
            raise FileNotFoundError(f"Missing package input: {source}")
        target = group_dir / source.name
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        records.append(
            {
                "group": group_name,
                "source_path": str(source),
                "package_path": str(target),
                "size_bytes": int(target.stat().st_size),
                "sha256": sha256_file(target),
            }
        )
    return records


def main():
    args = parse_args()
    package_dir = Path(args.package_dir)
    package_dir.mkdir(parents=True, exist_ok=True)

    records = []
    records.extend(copy_group(args.reference_files, package_dir, "reference"))
    records.extend(copy_group(args.asset_files, package_dir, "assets"))
    if args.metadata_files:
        records.extend(copy_group(args.metadata_files, package_dir, "metadata"))

    inventory_df = pd.DataFrame(records)
    inventory_path = Path(args.inventory_output)
    inventory_path.parent.mkdir(parents=True, exist_ok=True)
    inventory_df.to_csv(inventory_path, sep="\t", index=False)

    manifest = {
        "package_name": args.package_name,
        "package_dir": str(package_dir),
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "file_count": int(len(records)),
        "groups": inventory_df["group"].value_counts().to_dict(),
        "files": records,
    }
    Path(args.manifest_output).write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    readme_lines = [
        f"# Frozen Reference Package: {args.package_name}",
        "",
        "## Contents",
        "",
        f"- package_dir: `{package_dir}`",
        f"- created_at: `{manifest['created_at']}`",
        f"- file_count: `{manifest['file_count']}`",
        "",
        "## Groups",
        "",
    ]
    for group_name, count in manifest["groups"].items():
        readme_lines.append(f"- {group_name}: {count} files")
    Path(args.readme_output).write_text("\n".join(readme_lines) + "\n", encoding="utf-8")
    Path(args.done_output).write_text(f"{manifest['package_name']}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
