#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from pathlib import Path


def load_yaml_file(path_value: str | Path) -> dict:
    import yaml

    path = Path(path_value)
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def write_yaml_file(path_value: str | Path, payload: dict) -> None:
    import yaml

    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=True)


def resolve_existing_path(path_value: str | Path, label: str, project_root: Path | None = None, expect_dir: bool = False) -> Path:
    if not str(path_value).strip():
        raise ValueError(f"--{label} is required")

    candidate = Path(path_value)
    candidates = [candidate]
    if project_root is not None and not candidate.is_absolute():
        candidates.insert(0, project_root / candidate)

    for current in candidates:
        if current.exists():
            if expect_dir and not current.is_dir():
                raise ValueError(f"--{label} must point to a directory: {current}")
            if not expect_dir and not current.is_file():
                raise ValueError(f"--{label} must point to a file: {current}")
            return current.resolve()

    raise FileNotFoundError(f"--{label} not found: {path_value}")
