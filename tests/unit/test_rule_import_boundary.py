from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class RuleImportBoundaryTest(unittest.TestCase):
    def test_rules_do_not_import_scripts_namespace(self):
        rule_files = sorted((REPO_ROOT / "rules").glob("*.smk"))
        offending = []
        for path in rule_files:
            for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
                stripped = line.strip()
                if stripped.startswith("from scripts.") or stripped.startswith("import scripts."):
                    offending.append(f"{path.relative_to(REPO_ROOT)}:{line_number}:{stripped}")
        self.assertEqual(offending, [])

    def test_snakefile_does_not_extend_scripts_sys_path(self):
        snakefile = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertNotIn('sys.path.insert(0, str(PIPELINE_ROOT / "scripts"))', snakefile)

    def test_snakefile_uses_dedicated_script_entrypoint_include(self):
        snakefile = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn('include: "rules/script_entrypoints.smk"', snakefile)

    def test_snakefile_uses_reference_layout_include(self):
        snakefile = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn('include: "rules/reference_layout.smk"', snakefile)

    def test_snakefile_uses_predict_layout_include(self):
        snakefile = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn('include: "rules/predict_layout.smk"', snakefile)

    def test_snakefile_uses_runtime_layout_include(self):
        snakefile = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn('include: "rules/runtime_layout.smk"', snakefile)

    def test_snakefile_uses_pipeline_modes_include(self):
        snakefile = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn('include: "rules/pipeline_modes.smk"', snakefile)

    def test_snakefile_uses_target_assembly_include(self):
        snakefile = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn('include: "rules/target_assembly.smk"', snakefile)


if __name__ == "__main__":
    unittest.main()
