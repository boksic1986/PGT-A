from __future__ import annotations

from pathlib import Path


def test_wisecondorx_predict_command_includes_blacklist_when_configured():
    repo_root = Path(__file__).resolve().parents[2]
    rule_text = (repo_root / "rules" / "predict_workflow.smk").read_text(encoding="utf-8")

    assert "blacklist=([CNV_POSTPROCESS_BLACKLIST_BED] if CNV_POSTPROCESS_BLACKLIST_BED else [])" in rule_text
    assert "blacklist=CNV_POSTPROCESS_BLACKLIST_BED" in rule_text
    assert "blacklist_args=\"--blacklist {params.blacklist}\"" in rule_text
    assert "$blacklist_args" in rule_text
    assert "blacklist not passed because" in rule_text


def test_predict_cnv_uses_mask_only_npz_input_contract():
    repo_root = Path(__file__).resolve().parents[2]
    layout_text = (repo_root / "rules" / "predict_layout.smk").read_text(encoding="utf-8")
    rule_text = (repo_root / "rules" / "predict_workflow.smk").read_text(encoding="utf-8")

    assert 'CNV_PREPROCESS_STRATEGY = str(CNV_PREPROCESS_CFG.get("strategy", "mask_only"))' in layout_text
    assert "CNV_MASKED_NPZ" in layout_text
    assert "CNV_PREDICT_INPUT_NPZ = CNV_MASKED_NPZ if CNV_PREPROCESS_STRATEGY == \"mask_only\" else CNV_NPZ" in layout_text
    assert "rule cnv_mask_npz_for_predict" in rule_text
    assert "pgta.reference.preprocess mask_npz" in rule_text
    assert "npz=CNV_PREDICT_INPUT_NPZ" in rule_text
