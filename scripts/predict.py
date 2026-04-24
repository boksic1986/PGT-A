#!/biosoftware/miniconda/envs/snakemake_env/bin/python
try:
    from scripts._compat_entry import run_dispatcher
except ImportError:
    from _compat_entry import run_dispatcher

if __name__ == "__main__":
    raise SystemExit(run_dispatcher(__file__))
