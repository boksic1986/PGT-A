from __future__ import annotations

import tempfile
import unittest
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pgta.reference.cohort import load_selected_samples, parse_decisions


class ReferenceCohortResequencingTest(unittest.TestCase):
    def test_promoted_reseq_replaces_original(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            summary = tmp / "baseline_qc_summary.tsv"
            summary.write_text(
                "\n".join(
                    [
                        "sample_id\tqc_decision",
                        "S1\tPASS",
                        "S1_R2\tPASS",
                        "S2\tWARN",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            manifest = tmp / "resequencing_manifest.tsv"
            manifest.write_text(
                "\n".join(
                    [
                        "sample_id\tsource_sample_id\tlibrary_id\tfastq_r1\tfastq_r2\tsex_group\tbatch_group\treseq_reason\treplacement_policy\tstatus\tdecision_reason",
                        "S1_R2\tS1\tLIB1\t/r1.fq.gz\t/r2.fq.gz\tXX\tB1\tlow_qc\treplace_original\tpromoted\tqc passed",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            selected = load_selected_samples(
                summary,
                parse_decisions("PASS,WARN"),
                resequencing_manifest=manifest,
            )

        self.assertEqual(selected, ["S1_R2", "S2"])

    def test_candidate_reseq_is_not_selected_until_promoted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            summary = tmp / "baseline_qc_summary.tsv"
            summary.write_text(
                "\n".join(
                    [
                        "sample_id\tqc_decision",
                        "S3_R2\tPASS",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            manifest = tmp / "resequencing_manifest.tsv"
            manifest.write_text(
                "\n".join(
                    [
                        "sample_id\tsource_sample_id\tlibrary_id\tfastq_r1\tfastq_r2\tsex_group\tbatch_group\treseq_reason\treplacement_policy\tstatus\tdecision_reason",
                        "S3_R2\tS3\tLIB2\t/r1.fq.gz\t/r2.fq.gz\tXY\tB1\tlow_qc\treplace_original\tcandidate\tawaiting qc",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "No samples matched"):
                load_selected_samples(
                    summary,
                    parse_decisions("PASS"),
                    resequencing_manifest=manifest,
                )


if __name__ == "__main__":
    unittest.main()
