from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
GITATTRIBUTES = REPO_ROOT / ".gitattributes"


def _gitattributes_lines():
    assert GITATTRIBUTES.exists(), ".gitattributes must define workflow line endings"
    return [
        line.strip()
        for line in GITATTRIBUTES.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def test_snakemake_entrypoints_are_forced_to_lf():
    lines = _gitattributes_lines()

    required = {
        "Snakefile text eol=lf",
        "*.smk text eol=lf",
    }
    assert required.issubset(set(lines))


def test_workflow_support_files_are_forced_to_lf():
    lines = _gitattributes_lines()

    required = {
        "*.py text eol=lf",
        "*.yaml text eol=lf",
        "*.yml text eol=lf",
    }
    assert required.issubset(set(lines))
