from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    import pandas as pd
except ModuleNotFoundError as exc:  # pragma: no cover
    pd = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None
    from pgta.predict.branch_a import assemble_a_branch_candidates


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class BranchACandidateAssemblyTest(unittest.TestCase):
    def test_merges_adjacent_same_state_candidates_and_ranks_by_abs_z(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "S1_aberrations.bed"
            pd.DataFrame(
                [
                    {"chr": "7", "start": 1, "end": 100, "ratio": 0.10, "zscore": 6.0, "type": "gain"},
                    {"chr": "7", "start": 101, "end": 200, "ratio": 0.20, "zscore": 12.0, "type": "gain"},
                    {"chr": "8", "start": 1, "end": 150, "ratio": -0.10, "zscore": -8.0, "type": "loss"},
                ]
            ).to_csv(bed_path, sep="\t", index=False)

            candidates, raw = assemble_a_branch_candidates(bed_path, "S1")

        self.assertEqual(len(raw), 3)
        self.assertEqual(len(candidates), 2)
        top = candidates.iloc[0]
        self.assertEqual(top["chrom"], "chr7")
        self.assertEqual((int(top["start"]), int(top["end"])), (1, 200))
        self.assertEqual(top["state"], "gain")
        self.assertEqual(float(top["a_zscore"]), 12.0)
        self.assertEqual(int(top["a_source_event_count"]), 2)
        self.assertEqual(top["a_support_level"], "strong")
        self.assertEqual(top["candidate_id"], "S1.A0001")

    def test_does_not_merge_opposite_state_candidates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "S2_aberrations.bed"
            pd.DataFrame(
                [
                    {"chr": "13", "start": 1, "end": 100, "ratio": 0.10, "zscore": 11.0, "type": "gain"},
                    {"chr": "13", "start": 101, "end": 200, "ratio": -0.10, "zscore": -12.0, "type": "loss"},
                ]
            ).to_csv(bed_path, sep="\t", index=False)

            candidates, _raw = assemble_a_branch_candidates(bed_path, "S2")

        self.assertEqual(len(candidates), 2)
        self.assertEqual(set(candidates["state"].tolist()), {"gain", "loss"})


if __name__ == "__main__":
    unittest.main()
