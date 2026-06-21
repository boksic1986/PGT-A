from __future__ import annotations

import pandas as pd

from pgta.predict.branch_b.lowres_evidence import annotate_lowres_evidence


def test_2mb_same_direction_supports_3_to_4mb_candidate():
    candidates = pd.DataFrame(
        [
            {
                "sample_id": "H1",
                "candidate_id": "H1_chr8_gain",
                "chrom": "chr8",
                "start": 10_000_000,
                "end": 13_500_000,
                "state": "gain",
                "a_abs_zscore": 12.0,
            }
        ]
    )
    lowres_2mb = pd.DataFrame(
        [
            {
                "sample_id": "H1",
                "chrom": "chr8",
                "start": 10_000_000,
                "end": 14_000_000,
                "state": "gain",
                "a_abs_zscore": 6.4,
            }
        ]
    )

    annotated, summary = annotate_lowres_evidence(
        candidates,
        sample_id="H1",
        lowres_2mb_events=lowres_2mb,
    )
    row = annotated.iloc[0]

    assert row["lowres_2mb_support_label"] == "LOWRES_SAME_DIRECTION_SUPPORT"
    assert row["lowres_2mb_same_direction"] == 1
    assert row["lowres_2mb_overlap_fraction"] == 1.0
    assert row["lowres_consensus_label"] == "LOWRES_2MB_SUPPORT_FOR_3_4MB_EVENT"
    assert summary["lowres_2mb_support_label_counts"] == {"LOWRES_SAME_DIRECTION_SUPPORT": 1}


def test_3mb_no_support_for_short_weak_candidate_is_not_informative():
    candidates = pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6_chr21_gain",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
            }
        ]
    )

    annotated, _ = annotate_lowres_evidence(
        candidates,
        sample_id="H6",
        lowres_3mb_events=pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"]),
    )
    row = annotated.iloc[0]

    assert row["lowres_3mb_support_label"] == "LOWRES_NOT_INFORMATIVE_SHORT_OR_BOUNDARY_EVENT"
    assert row["lowres_3mb_same_direction"] == 0
    assert row["lowres_consensus_label"] == "LOWRES_CONTEXT_ONLY_SHORT_OR_BOUNDARY_EVENT"


def test_lowres_no_support_for_broad_candidate_is_context_not_filter():
    candidates = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": "S1_chr3_loss",
                "chrom": "chr3",
                "start": 1_000_000,
                "end": 9_500_000,
                "state": "loss",
                "a_abs_zscore": 8.0,
            }
        ]
    )

    annotated, _ = annotate_lowres_evidence(
        candidates,
        sample_id="S1",
        lowres_2mb_events=pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"]),
        lowres_3mb_events=pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"]),
    )
    row = annotated.iloc[0]

    assert row["lowres_2mb_support_label"] == "LOWRES_NO_SAME_DIRECTION_SUPPORT"
    assert row["lowres_3mb_support_label"] == "LOWRES_NO_SAME_DIRECTION_SUPPORT"
    assert row["lowres_consensus_label"] == "LOWRES_NO_SUPPORT_INFORMATIVE_BUT_NOT_FILTER"


def test_2mb_no_support_without_3mb_is_still_non_filter_context():
    candidates = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": "S1_chr5_gain",
                "chrom": "chr5",
                "start": 1_000_000,
                "end": 7_500_000,
                "state": "gain",
                "a_abs_zscore": 9.0,
            }
        ]
    )

    annotated, summary = annotate_lowres_evidence(
        candidates,
        sample_id="S1",
        lowres_2mb_events=pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"]),
    )
    row = annotated.iloc[0]

    assert row["lowres_2mb_support_label"] == "LOWRES_NO_SAME_DIRECTION_SUPPORT"
    assert row["lowres_3mb_support_label"] == "LOWRES_NOT_CONFIGURED"
    assert row["lowres_consensus_label"] == "LOWRES_NO_SUPPORT_INFORMATIVE_BUT_NOT_FILTER"
    assert summary["lowres_absence_is_filter_evidence"] is False
