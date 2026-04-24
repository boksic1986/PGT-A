from __future__ import annotations


BASELINE_QC_REQUIRED_COLUMNS = [
    "sample_id",
    "qc_decision",
    "qc_reason",
]

FRACTION_TRUTH_REQUIRED_COLUMNS = [
    "sample_id",
    "chrom",
    "expected_state",
    "abnormal_cell_fraction",
]

CNV_EVENT_KEY_COLUMNS = [
    "sample_id",
    "chrom",
    "start",
    "end",
    "state",
]
