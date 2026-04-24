import csv
import json
import math
import sqlite3
from datetime import datetime
from pathlib import Path


SUCCESS_MARKERS = (
    "completed successfully",
    "done",
    "finished",
    "written",
    "building gender reference",
    "common reference binsize written",
)
FAILURE_MARKERS = (
    "traceback",
    "error:",
    "exception",
    "workflowerror",
    "killed",
)


def ensure_parent(path_value):
    Path(path_value).parent.mkdir(parents=True, exist_ok=True)


def read_benchmark_row(path_value):
    with open(path_value, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty benchmark file: {path_value}")
    return rows[-1]


def parse_float(value):
    if value in (None, "", "NA", "nan"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_yaml_like(path_value):
    yaml = __import__("yaml")
    with open(path_value, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def detect_status(output_paths, log_path):
    outputs_exist = bool(output_paths) and all(Path(path).exists() for path in output_paths)
    if outputs_exist:
        return "success"
    if log_path and Path(log_path).exists():
        log_text = Path(log_path).read_text(encoding="utf-8", errors="ignore").lower()
        if any(marker in log_text for marker in FAILURE_MARKERS):
            return "log_failed"
        if any(marker in log_text for marker in SUCCESS_MARKERS):
            return "log_complete"
        return "benchmark_only"
    return "output_missing"


def enrich_parameters(record):
    binsize = record.get("binsize")
    pca_components = record.get("pca_components")
    parameter_path = record.get("parameter_path")
    if not parameter_path or not Path(parameter_path).exists():
        return binsize, pca_components

    suffixes = Path(parameter_path).suffixes
    if suffixes and suffixes[-1] in {".yaml", ".yml"}:
        payload = load_yaml_like(parameter_path)
        binsize = binsize if binsize is not None else payload.get("best_binsize") or payload.get("binsize")
        pca_components = (
            pca_components
            if pca_components is not None
            else payload.get("best_pca_components") or payload.get("pca_components")
        )
    elif Path(parameter_path).suffix == ".txt":
        text = Path(parameter_path).read_text(encoding="utf-8").strip()
        if text:
            try:
                binsize = binsize if binsize is not None else int(text)
            except ValueError:
                pass
    return binsize, pca_components


def percentile(values, q):
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[int(position)]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def serialize_output_path(output_paths):
    if not output_paths:
        return None
    if len(output_paths) == 1:
        return output_paths[0]
    return json.dumps(output_paths, ensure_ascii=True)


def create_tables(conn):
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS task_runs (
            id INTEGER PRIMARY KEY,
            module TEXT NOT NULL,
            rule_name TEXT NOT NULL,
            snakemake_target TEXT,
            sample_id TEXT,
            cohort_type TEXT,
            sample_bucket TEXT,
            binsize INTEGER,
            pca_components INTEGER,
            threads INTEGER,
            benchmark_path TEXT NOT NULL UNIQUE,
            log_path TEXT,
            output_path TEXT,
            status TEXT NOT NULL,
            wall_seconds REAL,
            max_rss_mb REAL,
            created_at TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS runtime_summary (
            module TEXT NOT NULL,
            rule_name TEXT NOT NULL,
            sample_bucket TEXT,
            threads INTEGER,
            binsize INTEGER,
            pca_components INTEGER,
            n_runs INTEGER NOT NULL,
            p50_minutes REAL,
            p80_minutes REAL,
            max_minutes REAL,
            recommended_check_after_minutes REAL,
            updated_at TEXT NOT NULL
        );
        """
    )


def collect_task_runs(records):
    task_rows = []
    for record in records:
        benchmark_path = record["benchmark_path"]
        if not Path(benchmark_path).exists():
            continue
        benchmark_row = read_benchmark_row(benchmark_path)
        wall_seconds = parse_float(benchmark_row.get("s"))
        max_rss_mb = parse_float(benchmark_row.get("max_rss"))
        binsize, pca_components = enrich_parameters(record)
        output_paths = record.get("output_paths", [])
        created_at = datetime.fromtimestamp(Path(benchmark_path).stat().st_mtime).isoformat(timespec="seconds")
        task_rows.append(
            {
                "module": record["module"],
                "rule_name": record["rule_name"],
                "snakemake_target": record.get("snakemake_target"),
                "sample_id": record.get("sample_id"),
                "cohort_type": record.get("cohort_type"),
                "sample_bucket": record.get("sample_bucket"),
                "binsize": binsize,
                "pca_components": pca_components,
                "threads": record.get("threads"),
                "benchmark_path": benchmark_path,
                "log_path": record.get("log_path"),
                "output_path": serialize_output_path(output_paths),
                "status": detect_status(output_paths, record.get("log_path")),
                "wall_seconds": wall_seconds,
                "max_rss_mb": max_rss_mb,
                "created_at": created_at,
            }
        )
    return task_rows


def write_task_runs(conn, task_rows):
    conn.executemany(
        """
        INSERT INTO task_runs (
            module, rule_name, snakemake_target, sample_id, cohort_type, sample_bucket,
            binsize, pca_components, threads, benchmark_path, log_path, output_path,
            status, wall_seconds, max_rss_mb, created_at
        ) VALUES (
            :module, :rule_name, :snakemake_target, :sample_id, :cohort_type, :sample_bucket,
            :binsize, :pca_components, :threads, :benchmark_path, :log_path, :output_path,
            :status, :wall_seconds, :max_rss_mb, :created_at
        )
        ON CONFLICT(benchmark_path) DO UPDATE SET
            module=excluded.module,
            rule_name=excluded.rule_name,
            snakemake_target=excluded.snakemake_target,
            sample_id=excluded.sample_id,
            cohort_type=excluded.cohort_type,
            sample_bucket=excluded.sample_bucket,
            binsize=excluded.binsize,
            pca_components=excluded.pca_components,
            threads=excluded.threads,
            log_path=excluded.log_path,
            output_path=excluded.output_path,
            status=excluded.status,
            wall_seconds=excluded.wall_seconds,
            max_rss_mb=excluded.max_rss_mb,
            created_at=excluded.created_at
        """,
        task_rows,
    )


def rebuild_runtime_summary(conn):
    rows = conn.execute(
        """
        SELECT module, rule_name, sample_bucket, threads, binsize, pca_components, wall_seconds
        FROM task_runs
        WHERE wall_seconds IS NOT NULL AND status IN ('success', 'log_complete')
        """
    ).fetchall()

    grouped = {}
    for row in rows:
        key = row[:6]
        grouped.setdefault(key, []).append(row[6] / 60.0)

    summary_rows = []
    updated_at = datetime.now().isoformat(timespec="seconds")
    for (module, rule_name, sample_bucket, threads, binsize, pca_components), minutes in grouped.items():
        p50 = percentile(minutes, 0.5)
        p80 = percentile(minutes, 0.8)
        max_minutes = max(minutes)
        recommended = None
        if p50 is not None and p80 is not None:
            recommended = math.ceil(max(p50 * 1.2, p80))
        summary_rows.append(
            {
                "module": module,
                "rule_name": rule_name,
                "sample_bucket": sample_bucket,
                "threads": threads,
                "binsize": binsize,
                "pca_components": pca_components,
                "n_runs": len(minutes),
                "p50_minutes": p50,
                "p80_minutes": p80,
                "max_minutes": max_minutes,
                "recommended_check_after_minutes": recommended,
                "updated_at": updated_at,
            }
        )

    conn.execute("DELETE FROM runtime_summary")
    conn.executemany(
        """
        INSERT INTO runtime_summary (
            module, rule_name, sample_bucket, threads, binsize, pca_components,
            n_runs, p50_minutes, p80_minutes, max_minutes, recommended_check_after_minutes, updated_at
        ) VALUES (
            :module, :rule_name, :sample_bucket, :threads, :binsize, :pca_components,
            :n_runs, :p50_minutes, :p80_minutes, :max_minutes, :recommended_check_after_minutes, :updated_at
        )
        """,
        summary_rows,
    )


def main():
    records = list(snakemake.params.records)
    db_path = snakemake.output.db
    done_path = snakemake.output.done

    ensure_parent(db_path)
    conn = sqlite3.connect(db_path)
    try:
        create_tables(conn)
        task_rows = collect_task_runs(records)
        write_task_runs(conn, task_rows)
        rebuild_runtime_summary(conn)
        conn.commit()
    finally:
        conn.close()

    Path(done_path).parent.mkdir(parents=True, exist_ok=True)
    Path(done_path).write_text(
        f"runtime tracking updated at {datetime.now().isoformat(timespec='seconds')}\n",
        encoding="utf-8",
    )

if __name__ == "__main__":
    main()
