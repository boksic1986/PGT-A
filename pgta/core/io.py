from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def ensure_parent(path_value: str | Path) -> Path:
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def read_tsv(path_value: str | Path, empty_ok: bool = False) -> pd.DataFrame:
    path = Path(path_value)
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    if df.empty and not empty_ok:
        raise ValueError(f"TSV is empty: {path}")
    return df


def write_tsv(path_value: str | Path, frame: pd.DataFrame) -> None:
    path = ensure_parent(path_value)
    frame.to_csv(path, sep="\t", index=False)


def read_json(path_value: str | Path, default=None):
    path = Path(path_value)
    if not path.exists():
        return {} if default is None else default
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path_value: str | Path, payload) -> None:
    path = ensure_parent(path_value)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
