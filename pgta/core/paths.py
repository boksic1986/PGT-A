from __future__ import annotations

from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def project_path(project_root: str | Path, *parts: str) -> Path:
    return Path(project_root).joinpath(*parts)


def resolve_project_path(project_root: str | Path, path_value: str | Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    return Path(project_root) / path
