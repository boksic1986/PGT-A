#!/biosoftware/miniconda/envs/snakemake_env/bin/python
try:
    from scripts._compat_entry import export_script
except ImportError:
    from _compat_entry import export_script

main = export_script(globals(), __file__)

if __name__ == "__main__" and main is not None:
    main()

