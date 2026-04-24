from __future__ import annotations

import importlib
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts import _compat_entry
from scripts._compat_entry import DISPATCHER_MODULE_MAP, SCRIPT_MODULE_MAP

OPTIONAL_LOCAL_DEPS = {"matplotlib", "numpy", "pandas"}
EXPECTED_SCRIPT_FILES = {
    "_compat_entry.py",
    "README.md",
    "init_pgta_project.py",
    "init_predict_config.py",
    "predict.py",
    "qc.py",
    "reference.py",
    "runtime.py",
    "validation.py",
}


class ScriptCompatTest(unittest.TestCase):
    def test_script_root_matches_expected_layout(self):
        script_files = {path.name for path in (REPO_ROOT / "scripts").iterdir() if path.is_file()}
        self.assertEqual(script_files, EXPECTED_SCRIPT_FILES)

    def test_wrapper_set_matches_map(self):
        wrappers = sorted(SCRIPT_MODULE_MAP)
        self.assertEqual(wrappers, sorted({"init_pgta_project.py", "init_predict_config.py", "runtime.py"}))

    def test_dispatcher_set_matches_map(self):
        dispatchers = sorted(DISPATCHER_MODULE_MAP)
        self.assertEqual(dispatchers, sorted({"predict.py", "qc.py", "reference.py", "validation.py"}))

    def test_python_entrypoints_use_clean_shebang_without_bom(self):
        for script_name in sorted(set(SCRIPT_MODULE_MAP) | set(DISPATCHER_MODULE_MAP)):
            with self.subTest(script=script_name):
                payload = (REPO_ROOT / "scripts" / script_name).read_bytes()
                self.assertFalse(payload.startswith(b"\xef\xbb\xbf"))
                self.assertTrue(payload.startswith(b"#!"))

    def test_wrappers_route_via_export_script(self):
        for script_name in SCRIPT_MODULE_MAP:
            with self.subTest(script=script_name):
                text = (REPO_ROOT / "scripts" / script_name).read_text(encoding="utf-8")
                self.assertIn("export_script(globals(), __file__)", text)

    def test_export_script_forwards_snakemake_into_target_module(self):
        fake_snakemake = object()
        fake_module = types.SimpleNamespace(
            __name__="pgta.core.runtime_tracking",
            __all__=["main"],
            main=lambda: 0,
        )

        with mock.patch.object(_compat_entry, "bootstrap_repo_root", return_value=REPO_ROOT):
            with mock.patch.object(_compat_entry.importlib, "import_module", return_value=fake_module):
                namespace = {"snakemake": fake_snakemake}
                main = _compat_entry.export_script(namespace, REPO_ROOT / "scripts" / "runtime.py")

        self.assertIs(main, fake_module.main)
        self.assertIs(getattr(fake_module, "snakemake"), fake_snakemake)

    def test_dispatchers_route_via_run_dispatcher(self):
        for script_name in DISPATCHER_MODULE_MAP:
            with self.subTest(script=script_name):
                text = (REPO_ROOT / "scripts" / script_name).read_text(encoding="utf-8")
                self.assertIn("run_dispatcher(__file__)", text)

    def test_dispatcher_bootstraps_repo_root_before_import(self):
        calls = []

        def fake_bootstrap():
            calls.append("bootstrap")
            return REPO_ROOT

        def fake_import_module(module_name):
            calls.append(("import", module_name))
            return types.SimpleNamespace(main=lambda: 0, __name__=module_name)

        with mock.patch.object(_compat_entry, "bootstrap_repo_root", side_effect=fake_bootstrap):
            with mock.patch.object(_compat_entry.importlib, "import_module", side_effect=fake_import_module):
                result = _compat_entry.run_dispatcher(REPO_ROOT / "scripts" / "predict.py", ["cnv_calling"])

        self.assertEqual(result, 0)
        self.assertGreaterEqual(len(calls), 2)
        self.assertEqual(calls[0], "bootstrap")
        self.assertEqual(calls[1], ("import", "pgta.predict.branch_b.calling"))

    def test_mapped_modules_import(self):
        all_mappings = dict(SCRIPT_MODULE_MAP)
        for action_map in DISPATCHER_MODULE_MAP.values():
            all_mappings.update(action_map)
        for entry_name, module_name in sorted(all_mappings.items()):
            with self.subTest(entry=entry_name, module=module_name):
                try:
                    importlib.import_module(module_name)
                except ModuleNotFoundError as exc:
                    if exc.name in OPTIONAL_LOCAL_DEPS:
                        self.skipTest(
                            f"local environment missing optional dependency: {exc.name}"
                        )
                    raise


if __name__ == "__main__":
    unittest.main()
