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
except ModuleNotFoundError as exc:  # pragma: no cover - environment-dependent skip
    pd = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None
    import pgta.predict.report as report_module
    from pgta.predict.report import (
        build_plot_lookup,
        format_biological_candidate_conclusion,
        format_technical_conclusion,
        load_branch_a_validation_summaries,
        load_branch_b_evidence_summaries,
        load_branch_s_summaries,
        summarize_branch_a_validation_gate,
        summarize_branch_b_events,
    )


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class CnvReportRankingTest(unittest.TestCase):
    def test_build_plot_lookup_indexes_existing_svg_by_sample_id(self):
        lookup = build_plot_lookup(
            [
                "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y1.final_cnv.svg",
                "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y2.branch_ab.svg",
            ]
        )

        self.assertEqual(
            lookup,
            {
                "Y1": "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y1.final_cnv.svg",
                "Y2": "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y2.branch_ab.svg",
            },
        )

    def test_sex_review_event_is_suppressed_when_autosomal_kept_event_exists(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "Y6",
                    "event_id": "event_sex",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 12.5,
                    "n_bins": 20,
                    "chrom": "chrX",
                    "start": 11_000_000,
                    "end": 37_000_000,
                    "state": "loss",
                    "technical_confidence": "moderate",
                    "artifact_flags": "sex_chromosome_event",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                },
                {
                    "sample_id": "Y6",
                    "event_id": "event_auto",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 9.2,
                    "n_bins": 12,
                    "chrom": "chr7",
                    "start": 63_000_001,
                    "end": 159_000_000,
                    "state": "gain",
                    "technical_confidence": "moderate",
                    "artifact_flags": "low_level_signal",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                },
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_suppressed_sex_review_events"]), 1)
        self.assertEqual(row["branch_b_top_event"], "chr7:63000001-159000000 gain [review/moderate]")

    def test_only_sex_review_event_remains_displayable(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "Y3",
                    "event_id": "event_sex",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 9.5,
                    "n_bins": 18,
                    "chrom": "chrX",
                    "start": 11_000_000,
                    "end": 37_000_000,
                    "state": "loss",
                    "technical_confidence": "moderate",
                    "artifact_flags": "sex_chromosome_event",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                }
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_suppressed_sex_review_events"]), 0)
        self.assertEqual(row["branch_b_top_event"], "chrX:11000000-37000000 loss [review/moderate]")

        conclusion_row = row.to_dict()
        conclusion_row["a_branch_top_event"] = ""
        technical = format_technical_conclusion(conclusion_row)
        biological = format_biological_candidate_conclusion(conclusion_row)

        self.assertIn("chrX:11000000-37000000 loss", technical)
        self.assertEqual(biological, "chrX:11000000-37000000 loss [review/moderate]")

    def test_pass_sex_chromosome_event_remains_displayable(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "S1",
                    "event_id": "event_sex_pass",
                    "keep_event": 1,
                    "artifact_status": "pass",
                    "priority_score": 15.0,
                    "n_bins": 30,
                    "chrom": "chrX",
                    "start": 5_000_000,
                    "end": 40_000_000,
                    "state": "loss",
                    "technical_confidence": "high",
                    "artifact_flags": "sex_chromosome_event",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "clean_support_and_statistical_support",
                }
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_suppressed_sex_review_events"]), 0)
        self.assertEqual(row["branch_b_top_event"], "chrX:5000000-40000000 loss [pass/high]")

    def test_cnvseq_reportable_tier_is_counted_and_prioritized_for_top_event(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "H1",
                    "event_id": "small_high_priority",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 100.0,
                    "n_bins": 2,
                    "chrom": "chr8",
                    "start": 10_000_000,
                    "end": 11_500_000,
                    "state": "loss",
                    "technical_confidence": "moderate",
                    "artifact_flags": "review_1_2mb",
                    "downgrade_reason": "cnvseq_subreportable_size_review",
                    "filter_reason": "",
                    "retain_reason": "",
                    "cnvseq_report_tier": "review_1_2mb",
                    "cnvseq_reportable": 0,
                },
                {
                    "sample_id": "H1",
                    "event_id": "large_lower_priority",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 10.0,
                    "n_bins": 4,
                    "chrom": "chr16",
                    "start": 46_000_001,
                    "end": 90_000_000,
                    "state": "gain",
                    "technical_confidence": "moderate",
                    "artifact_flags": "",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                    "cnvseq_report_tier": "reportable",
                    "cnvseq_reportable": 1,
                },
                {
                    "sample_id": "H1",
                    "event_id": "artifact",
                    "keep_event": 0,
                    "artifact_status": "artifact",
                    "priority_score": 999.0,
                    "n_bins": 20,
                    "chrom": "chr19",
                    "start": 1,
                    "end": 3_000_000,
                    "state": "gain",
                    "technical_confidence": "rejected",
                    "artifact_flags": "",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                    "cnvseq_report_tier": "reportable",
                    "cnvseq_reportable": 1,
                },
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_reportable_events"]), 1)
        self.assertEqual(int(row["branch_b_review_tier_events"]), 1)
        self.assertEqual(row["branch_b_top_event"], "chr16:46000001-90000000 gain [review/moderate; reportable]")

    def test_report_loads_p3_and_p5_summary_only_status(self):
        with tempfile.TemporaryDirectory() as temp_dir_value:
            temp_dir = Path(temp_dir_value)
            branch_b_summary = temp_dir / "S1.branch_b.summary.json"
            branch_b_summary.write_text(
                """
                {
                  "sample_id": "S1",
                  "candidate_count": 4,
                  "p3_disposition_counts": {"REVIEW_REQUIRED": 3, "NO_CALL_CONTRACT_RISK": 1},
                  "review_burden_count": 3,
                  "missing_evidence_candidate_count": 1,
                  "background_status_counts": {"UNKNOWN_BACKGROUND": 4},
                  "final_report_impact": "none_shadow_only"
                }
                """,
                encoding="utf-8",
            )
            branch_s_summary = temp_dir / "S1.branch_s.summary.json"
            branch_s_summary.write_text(
                """
                {
                  "sample_id": "S1",
                  "sca_candidate_state": "X_LOSS",
                  "sca_confidence_tier": "SCA_REVIEW_WEAK",
                  "sca_report_layer_class": "sca_internal_review_event",
                  "sca_report_layer_reason": "sca_review_weak_with_nonpar_corroboration",
                  "sca_output_mode": "review_development_only",
                  "report_text_status": "development_only_not_final_reportable",
                  "sca_uncertainty_reason": "locked_sca_truth_incomplete"
                }
                """,
                encoding="utf-8",
            )

            branch_b_df = load_branch_b_evidence_summaries([str(branch_b_summary)])
            branch_s_df = load_branch_s_summaries([str(branch_s_summary)])

        self.assertEqual(int(branch_b_df.iloc[0]["branch_b_evidence_candidate_count"]), 4)
        self.assertEqual(int(branch_b_df.iloc[0]["branch_b_evidence_review_required_count"]), 3)
        self.assertEqual(branch_b_df.iloc[0]["branch_b_evidence_background_status"], "UNKNOWN_BACKGROUND=4")
        self.assertEqual(branch_b_df.iloc[0]["branch_b_evidence_final_report_impact"], "none_shadow_only")
        self.assertEqual(branch_s_df.iloc[0]["sca_candidate_state"], "X_LOSS")
        self.assertEqual(branch_s_df.iloc[0]["sca_output_mode"], "review_development_only")
        self.assertEqual(branch_s_df.iloc[0]["sca_report_layer_class"], "sca_internal_review_event")
        self.assertEqual(
            branch_s_df.iloc[0]["sca_report_layer_reason"],
            "sca_review_weak_with_nonpar_corroboration",
        )

    def test_report_loads_branch_b_v2_burden_summary_for_development_display(self):
        self.assertTrue(hasattr(report_module, "load_branch_b_v2_burden_summaries"))

        with tempfile.TemporaryDirectory() as temp_dir_value:
            temp_dir = Path(temp_dir_value)
            sample_summary = temp_dir / "sample_summary.tsv"
            sample_summary.write_text(
                "\t".join(
                    [
                        "sample_id",
                        "candidate_count",
                        "v2_report_event_count",
                        "v2_internal_review_event_count",
                        "v2_filtered_event_count",
                        "v2_branch_s_event_count",
                        "v2_report_candidate_burden_count",
                        "v2_review_candidate_burden_count",
                        "v2_background_unknown_review_burden_count",
                        "v2_technical_risk_burden_count",
                        "v2_branch_s_review_burden_count",
                        "sample_report_burden_flag",
                        "sample_report_burden_reason",
                    ]
                )
                + "\n"
                + "\t".join(["S1", "5", "1", "2", "1", "1", "1", "3", "2", "1", "1", "1", "report_event_count_ge_3"])
                + "\n",
                encoding="utf-8",
            )
            benchmark_summary = temp_dir / "summary.json"
            benchmark_summary.write_text(
                """
                {
                  "status": "ready",
                  "reference_id": "h_r0_shadow_ref_20260619",
                  "sample_count": 1,
                  "candidate_count": 5,
                  "v2_report_event_count": 1,
                  "v2_internal_review_event_count": 2,
                  "v2_filtered_event_count": 1,
                  "v2_branch_s_event_count": 1,
                  "v2_report_candidate_burden_count": 1,
                  "v2_review_candidate_burden_count": 3,
                  "v2_background_unknown_review_burden_count": 2,
                  "v2_technical_risk_burden_count": 1,
                  "v2_branch_s_review_burden_count": 1,
                  "sample_report_burden_flag_count": 1,
                  "sample_report_burden_threshold": 3,
                  "legacy_branch_b_decision_fields_used": false,
                  "final_report_impact": "none_shadow_only",
                  "benchmark_scope": "v2_classifier_rows_only"
                }
                """,
                encoding="utf-8",
            )

            v2_df, v2_contract = report_module.load_branch_b_v2_burden_summaries(
                sample_summary_path=str(sample_summary),
                benchmark_summary_path=str(benchmark_summary),
                report_reference_id="h_r0_shadow_ref_20260619",
            )

        self.assertEqual(v2_contract["branch_b_v2_burden_status"], "ready")
        self.assertEqual(v2_contract["branch_b_v2_legacy_fields_used"], False)
        self.assertEqual(v2_contract["branch_b_v2_final_impact"], "development_review_only")
        self.assertEqual(v2_contract["branch_b_v2_report_event_count"], 1)
        self.assertEqual(v2_contract["branch_b_v2_internal_review_event_count"], 2)
        self.assertEqual(v2_contract["branch_b_v2_filtered_event_count"], 1)
        self.assertEqual(v2_contract["branch_b_v2_branch_s_event_count"], 1)
        self.assertEqual(v2_contract["branch_b_v2_sample_report_burden_flag_count"], 1)
        self.assertEqual(v2_contract["branch_b_v2_sample_report_burden_threshold"], 3)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_report_candidate_count"]), 1)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_report_event_count"]), 1)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_internal_review_event_count"]), 2)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_filtered_event_count"]), 1)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_branch_s_event_count"]), 1)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_background_unknown_review_count"]), 2)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_branch_s_review_count"]), 1)
        self.assertEqual(int(v2_df.iloc[0]["branch_b_v2_sample_report_burden_flag"]), 1)
        self.assertEqual(v2_df.iloc[0]["branch_b_v2_sample_report_burden_reason"], "report_event_count_ge_3")

    def test_technical_conclusion_marks_branch_b_v2_development_only_burden(self):
        row = pd.Series(
            {
                "branch_b_kept_events": 0,
                "branch_b_suppressed_sex_review_events": 0,
                "a_branch_event_count": 1,
                "a_branch_top_event": "chr21:20700000-22300000 gain",
                "a_branch_review_candidate_count": 1,
                "a_branch_review_shortlist": "chr21:20700000-22300000 gain z=7.11",
                "a_branch_strong_signal_count": 0,
                "branch_b_v2_burden_status": "ready",
                "branch_b_v2_background_unknown_review_count": 2,
                "branch_b_v2_branch_s_review_count": 1,
                "branch_b_v2_technical_risk_review_count": 0,
                "branch_b_v2_report_candidate_count": 0,
                "branch_b_v2_report_event_count": 1,
                "branch_b_v2_internal_review_event_count": 2,
                "branch_b_v2_filtered_event_count": 3,
                "branch_b_v2_branch_s_event_count": 1,
                "branch_b_v2_final_impact": "development_review_only",
                "branch_b_v2_legacy_fields_used": False,
                "qc_status": "PASS",
                "sex_call": "XY",
            }
        )

        conclusion = format_technical_conclusion(row)

        self.assertIn("Legacy/current-code Branch B kept 0 events", conclusion)
        self.assertIn("Bv2_burden_status=ready", conclusion)
        self.assertIn("Bv2_background_unknown_review=2", conclusion)
        self.assertIn("Bv2_branch_s_review=1", conclusion)
        self.assertIn("Bv2_report_candidate=0", conclusion)
        self.assertIn("Bv2_report_event=1", conclusion)
        self.assertIn("Bv2_internal_review_event=2", conclusion)
        self.assertIn("Bv2_filtered_event=3", conclusion)
        self.assertIn("Bv2_branch_s_event=1", conclusion)
        self.assertIn("Bv2_final_impact=development_review_only", conclusion)
        self.assertIn("Bv2_legacy_fields_used=False", conclusion)

    def test_biological_conclusion_prefers_v2_report_layer_over_legacy_branch_b_top_event(self):
        row = pd.Series(
            {
                "branch_b_top_event": "legacy chr16:1-5000000 gain [review/high]",
                "branch_b_v2_burden_status": "ready",
                "branch_b_v2_report_event_count": 2,
                "branch_b_v2_internal_review_event_count": 3,
                "branch_b_v2_filtered_event_count": 4,
                "branch_b_v2_branch_s_event_count": 1,
                "branch_b_v2_final_impact": "development_review_only",
            }
        )

        conclusion = format_biological_candidate_conclusion(row)

        self.assertIn("Branch B V2 report-layer", conclusion)
        self.assertIn("report_events=2", conclusion)
        self.assertIn("internal_review_events=3", conclusion)
        self.assertIn("filtered_events=4", conclusion)
        self.assertIn("branch_s_events=1", conclusion)
        self.assertNotIn("legacy chr16", conclusion)

    def test_technical_conclusion_marks_p6_review_development_context(self):
        row = pd.Series(
            {
                "branch_b_kept_events": 0,
                "branch_b_suppressed_sex_review_events": 0,
                "a_branch_event_count": 2,
                "a_branch_top_event": "chr21:20700000-22300000 gain",
                "a_branch_review_candidate_count": 2,
                "a_branch_review_shortlist": "chr21:20700000-22300000 gain z=7.11",
                "a_branch_strong_signal_count": 0,
                "branch_b_evidence_candidate_count": 2,
                "branch_b_evidence_review_required_count": 2,
                "branch_b_evidence_background_status": "UNKNOWN_BACKGROUND=2",
                "branch_b_evidence_final_report_impact": "none_shadow_only",
                "sca_candidate_state": "X_LOSS",
                "sca_confidence_tier": "SCA_REVIEW_WEAK",
                "sca_report_layer_class": "sca_internal_review_event",
                "sca_report_layer_reason": "sca_review_weak_with_nonpar_corroboration",
                "sca_output_mode": "review_development_only",
                "report_text_status": "development_only_not_final_reportable",
                "qc_status": "PASS",
                "sex_call": "XX",
            }
        )

        conclusion = format_technical_conclusion(row)

        self.assertIn("P3_branch_b_evidence_candidates=2", conclusion)
        self.assertIn("P3_review_required=2", conclusion)
        self.assertIn("P3_background=UNKNOWN_BACKGROUND=2", conclusion)
        self.assertIn("P3_report_impact=none_shadow_only", conclusion)
        self.assertIn("SCA_mode=review_development_only", conclusion)
        self.assertIn("SCA_report_layer=sca_internal_review_event", conclusion)
        self.assertIn("SCA_status=development_only_not_final_reportable", conclusion)

    def test_report_summarizes_branch_a_no_fn_gate_for_same_reference(self):
        with tempfile.TemporaryDirectory() as temp_dir_value:
            temp_dir = Path(temp_dir_value)
            y_summary = temp_dir / "y.branch_a.summary.json"
            y_summary.write_text(
                """
                {
                  "reference_id": "h_r0_shadow_ref_20260619",
                  "sample_count": 8,
                  "truth_event_count": 10,
                  "truth_detected_count": 10,
                  "FN_count": 0,
                  "H6_chr21_status": "not_applicable"
                }
                """,
                encoding="utf-8",
            )
            h_summary = temp_dir / "h.branch_a.summary.json"
            h_summary.write_text(
                """
                {
                  "reference_id": "h_r0_shadow_ref_20260619",
                  "sample_count": 16,
                  "truth_event_count": 10,
                  "truth_detected_count": 10,
                  "FN_count": 0,
                  "H6_chr21_status": "detected"
                }
                """,
                encoding="utf-8",
            )

            summaries = load_branch_a_validation_summaries([str(y_summary), str(h_summary)])
            gate = summarize_branch_a_validation_gate("h_r0_shadow_ref_20260619", summaries)

        self.assertEqual(gate["status"], "passed_no_fn")
        self.assertEqual(gate["summary_count"], 2)
        self.assertEqual(gate["truth_event_count"], 20)
        self.assertEqual(gate["truth_detected_count"], 20)
        self.assertEqual(gate["FN_count"], 0)
        self.assertEqual(gate["H6_chr21_status"], "detected")
        self.assertEqual(gate["same_reference_config_status"], "matched")

    def test_report_flags_branch_a_validation_reference_mismatch(self):
        with tempfile.TemporaryDirectory() as temp_dir_value:
            temp_dir = Path(temp_dir_value)
            summary_path = temp_dir / "old.branch_a.summary.json"
            summary_path.write_text(
                """
                {
                  "reference_id": "old_ref",
                  "truth_event_count": 10,
                  "truth_detected_count": 10,
                  "FN_count": 0
                }
                """,
                encoding="utf-8",
            )

            summaries = load_branch_a_validation_summaries([str(summary_path)])
            gate = summarize_branch_a_validation_gate("h_r0_shadow_ref_20260619", summaries)

        self.assertEqual(gate["status"], "reference_mismatch")
        self.assertEqual(gate["same_reference_config_status"], "mismatch")


if __name__ == "__main__":
    unittest.main()
