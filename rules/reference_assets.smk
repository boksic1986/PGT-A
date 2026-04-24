if REFERENCE_ASSET_TARGET_FILES:
    if TUNING_ENABLED and "reference_qc" in REQUESTED_TARGETS:
        rule build_reference_bin_annotations:
            input:
                best=TUNING_BEST,
                metadata=RUN_METADATA
            output:
                atomic_bins=REFERENCE_ATOMIC_BINS,
                analysis_bins=REFERENCE_ANALYSIS_BINS,
                qc_bins=REFERENCE_QC_BINS,
                atomic_annotations=REFERENCE_ATOMIC_BIN_ANNOTATIONS,
                analysis_annotations=REFERENCE_ANALYSIS_BIN_ANNOTATIONS,
                qc_annotations=REFERENCE_QC_BIN_ANNOTATIONS,
                summary_json=REFERENCE_BIN_SUMMARY_JSON
            log:
                project_path("logs", "reference_assets", "build_bin_annotations.log")
            params:
                python_bin=config["biosoft"]["python"],
                script=SCRIPT_BUILD_BIN_ANNOTATIONS,
                script_action=SCRIPT_BUILD_BIN_ANNOTATIONS_ACTION,
                reference_fasta=config["core"]["reference_genome"],
                chromosomes=" ".join(config["core"]["chromosome_list"]),
                atomic_binsize=ATOMIC_BINSIZE,
                qc_binsize=QC_BINSIZE,
                analysis_binsize=lambda wildcards, input: read_best_binsize_from_yaml(input.best, ANALYSIS_BINSIZE_DEFAULT),
                par_regions=" ".join(CNV_POSTPROCESS_PAR_REGIONS),
                xtr_regions=" ".join(CNV_POSTPROCESS_XTR_REGIONS),
                segmental_duplication_bed=CNV_POSTPROCESS_SEGMENTAL_DUPLICATION_BED,
                low_mappability_bed=CNV_POSTPROCESS_LOW_MAPPABILITY_BED,
                gap_bed=CNV_POSTPROCESS_GAP_BED,
                repeat_bed=CNV_POSTPROCESS_REPEAT_BED,
                blacklist_bed=CNV_POSTPROCESS_BLACKLIST_BED,
                sex_homology_bed=CNV_POSTPROCESS_SEX_HOMOLOGY_BED,
                ambiguous_alignment_bed=CNV_POSTPROCESS_AMBIGUOUS_ALIGNMENT_BED
            threads: 1
            shell:
                r"""
                mkdir -p "$(dirname {output.atomic_bins})" "$(dirname {log})"
                cmd=(
                    {params.python_bin:q} {params.script:q} {params.script_action:q} annotations \
                    --reference-fasta {params.reference_fasta:q} \
                    --chromosomes {params.chromosomes} \
                    --atomic-bin-size {params.atomic_binsize} \
                    --analysis-bin-size {params.analysis_binsize} \
                    --qc-bin-size {params.qc_binsize} \
                    --atomic-bins-output {output.atomic_bins:q} \
                    --analysis-bins-output {output.analysis_bins:q} \
                    --qc-bins-output {output.qc_bins:q} \
                    --atomic-annotations-output {output.atomic_annotations:q} \
                    --analysis-annotations-output {output.analysis_annotations:q} \
                    --qc-annotations-output {output.qc_annotations:q} \
                    --summary-json-output {output.summary_json:q} \
                    --log {log:q}
                )
                if [ -n "{params.par_regions}" ]; then
                    for region in {params.par_regions}; do
                        cmd+=(--par-region "$region")
                    done
                fi
                if [ -n "{params.xtr_regions}" ]; then
                    for region in {params.xtr_regions}; do
                        cmd+=(--xtr-region "$region")
                    done
                fi
                if [ -n "{params.segmental_duplication_bed}" ]; then cmd+=(--segmental-duplication-bed {params.segmental_duplication_bed:q}); fi
                if [ -n "{params.low_mappability_bed}" ]; then cmd+=(--low-mappability-bed {params.low_mappability_bed:q}); fi
                if [ -n "{params.gap_bed}" ]; then cmd+=(--gap-centromere-telomere-bed {params.gap_bed:q}); fi
                if [ -n "{params.repeat_bed}" ]; then cmd+=(--repeat-rich-bed {params.repeat_bed:q}); fi
                if [ -n "{params.blacklist_bed}" ]; then cmd+=(--blacklist-bed {params.blacklist_bed:q}); fi
                if [ -n "{params.sex_homology_bed}" ]; then cmd+=(--sex-homology-bed {params.sex_homology_bed:q}); fi
                if [ -n "{params.ambiguous_alignment_bed}" ]; then cmd+=(--ambiguous-alignment-bed {params.ambiguous_alignment_bed:q}); fi
                "${{cmd[@]}}" >> {log:q} 2>&1
                """
    else:
        rule build_reference_bin_annotations:
            input:
                metadata=RUN_METADATA
            output:
                atomic_bins=REFERENCE_ATOMIC_BINS,
                analysis_bins=REFERENCE_ANALYSIS_BINS,
                qc_bins=REFERENCE_QC_BINS,
                atomic_annotations=REFERENCE_ATOMIC_BIN_ANNOTATIONS,
                analysis_annotations=REFERENCE_ANALYSIS_BIN_ANNOTATIONS,
                qc_annotations=REFERENCE_QC_BIN_ANNOTATIONS,
                summary_json=REFERENCE_BIN_SUMMARY_JSON
            log:
                project_path("logs", "reference_assets", "build_bin_annotations.log")
            params:
                python_bin=config["biosoft"]["python"],
                script=SCRIPT_BUILD_BIN_ANNOTATIONS,
                script_action=SCRIPT_BUILD_BIN_ANNOTATIONS_ACTION,
                reference_fasta=config["core"]["reference_genome"],
                chromosomes=" ".join(config["core"]["chromosome_list"]),
                atomic_binsize=ATOMIC_BINSIZE,
                analysis_binsize=ANALYSIS_BINSIZE_DEFAULT,
                qc_binsize=QC_BINSIZE,
                par_regions=" ".join(CNV_POSTPROCESS_PAR_REGIONS),
                xtr_regions=" ".join(CNV_POSTPROCESS_XTR_REGIONS),
                segmental_duplication_bed=CNV_POSTPROCESS_SEGMENTAL_DUPLICATION_BED,
                low_mappability_bed=CNV_POSTPROCESS_LOW_MAPPABILITY_BED,
                gap_bed=CNV_POSTPROCESS_GAP_BED,
                repeat_bed=CNV_POSTPROCESS_REPEAT_BED,
                blacklist_bed=CNV_POSTPROCESS_BLACKLIST_BED,
                sex_homology_bed=CNV_POSTPROCESS_SEX_HOMOLOGY_BED,
                ambiguous_alignment_bed=CNV_POSTPROCESS_AMBIGUOUS_ALIGNMENT_BED
            threads: 1
            shell:
                r"""
                mkdir -p "$(dirname {output.atomic_bins})" "$(dirname {log})"
                cmd=(
                    {params.python_bin:q} {params.script:q} {params.script_action:q} annotations \
                    --reference-fasta {params.reference_fasta:q} \
                    --chromosomes {params.chromosomes} \
                    --atomic-bin-size {params.atomic_binsize} \
                    --analysis-bin-size {params.analysis_binsize} \
                    --qc-bin-size {params.qc_binsize} \
                    --atomic-bins-output {output.atomic_bins:q} \
                    --analysis-bins-output {output.analysis_bins:q} \
                    --qc-bins-output {output.qc_bins:q} \
                    --atomic-annotations-output {output.atomic_annotations:q} \
                    --analysis-annotations-output {output.analysis_annotations:q} \
                    --qc-annotations-output {output.qc_annotations:q} \
                    --summary-json-output {output.summary_json:q} \
                    --log {log:q}
                )
                if [ -n "{params.par_regions}" ]; then
                    for region in {params.par_regions}; do
                        cmd+=(--par-region "$region")
                    done
                fi
                if [ -n "{params.xtr_regions}" ]; then
                    for region in {params.xtr_regions}; do
                        cmd+=(--xtr-region "$region")
                    done
                fi
                if [ -n "{params.segmental_duplication_bed}" ]; then cmd+=(--segmental-duplication-bed {params.segmental_duplication_bed:q}); fi
                if [ -n "{params.low_mappability_bed}" ]; then cmd+=(--low-mappability-bed {params.low_mappability_bed:q}); fi
                if [ -n "{params.gap_bed}" ]; then cmd+=(--gap-centromere-telomere-bed {params.gap_bed:q}); fi
                if [ -n "{params.repeat_bed}" ]; then cmd+=(--repeat-rich-bed {params.repeat_bed:q}); fi
                if [ -n "{params.blacklist_bed}" ]; then cmd+=(--blacklist-bed {params.blacklist_bed:q}); fi
                if [ -n "{params.sex_homology_bed}" ]; then cmd+=(--sex-homology-bed {params.sex_homology_bed:q}); fi
                if [ -n "{params.ambiguous_alignment_bed}" ]; then cmd+=(--ambiguous-alignment-bed {params.ambiguous_alignment_bed:q}); fi
                "${{cmd[@]}}" >> {log:q} 2>&1
                """

    if "baseline_qc" in REQUESTED_TARGETS:
        rule build_reference_masks:
            input:
                atomic_annotations=REFERENCE_ATOMIC_BIN_ANNOTATIONS,
                analysis_annotations=REFERENCE_ANALYSIS_BIN_ANNOTATIONS,
                qc_annotations=REFERENCE_QC_BIN_ANNOTATIONS,
                profiles=expand(BASELINE_QC_PROFILE_TSV, sample=BASELINE_SAMPLE_IDS),
                metadata=RUN_METADATA
            output:
                hard_mask=REFERENCE_HARD_MASK_TSV,
                soft_mask=REFERENCE_SOFT_MASK_TSV,
                dynamic_mask=REFERENCE_DYNAMIC_MASK_TSV,
                combined_mask=REFERENCE_COMBINED_MASK_TSV,
                hard_json=REFERENCE_HARD_MASK_JSON,
                soft_json=REFERENCE_SOFT_MASK_JSON,
                dynamic_json=REFERENCE_DYNAMIC_MASK_JSON,
                combined_json=REFERENCE_COMBINED_MASK_JSON,
                summary_json=REFERENCE_MASK_SUMMARY_JSON
            log:
                project_path("logs", "reference_assets", "build_bin_mask.log")
            params:
                python_bin=config["biosoft"]["python"],
                script=SCRIPT_BUILD_BIN_ANNOTATIONS,
                script_action=SCRIPT_BUILD_BIN_ANNOTATIONS_ACTION,
                hard_n_fraction=HARD_MASK_N_FRACTION,
                soft_gc_low=SOFT_MASK_GC_LOW,
                soft_gc_high=SOFT_MASK_GC_HIGH,
                dynamic_z_frac=DYNAMIC_MASK_Z_FRAC,
                dynamic_median_abs_z=DYNAMIC_MASK_MEDIAN_ABS_Z
            threads: 1
            shell:
                r"""
                mkdir -p "$(dirname {output.hard_mask})" "$(dirname {log})"
                {params.python_bin:q} {params.script:q} {params.script_action:q} \
                    mask \
                    --annotation-tsvs {input.atomic_annotations:q} {input.analysis_annotations:q} {input.qc_annotations:q} \
                    --profile-tsvs {input.profiles:q} \
                    --hard-mask-output {output.hard_mask:q} \
                    --soft-mask-output {output.soft_mask:q} \
                    --dynamic-mask-output {output.dynamic_mask:q} \
                    --combined-mask-output {output.combined_mask:q} \
                    --hard-mask-json-output {output.hard_json:q} \
                    --soft-mask-json-output {output.soft_json:q} \
                    --dynamic-mask-json-output {output.dynamic_json:q} \
                    --combined-mask-json-output {output.combined_json:q} \
                    --summary-json-output {output.summary_json:q} \
                    --hard-n-fraction {params.hard_n_fraction} \
                    --soft-gc-low {params.soft_gc_low} \
                    --soft-gc-high {params.soft_gc_high} \
                    --dynamic-z-frac-threshold {params.dynamic_z_frac} \
                    --dynamic-median-abs-z-threshold {params.dynamic_median_abs_z} \
                    >> {log:q} 2>&1
                """
    else:
        rule build_reference_masks:
            input:
                atomic_annotations=REFERENCE_ATOMIC_BIN_ANNOTATIONS,
                analysis_annotations=REFERENCE_ANALYSIS_BIN_ANNOTATIONS,
                qc_annotations=REFERENCE_QC_BIN_ANNOTATIONS,
                metadata=RUN_METADATA
            output:
                hard_mask=REFERENCE_HARD_MASK_TSV,
                soft_mask=REFERENCE_SOFT_MASK_TSV,
                dynamic_mask=REFERENCE_DYNAMIC_MASK_TSV,
                combined_mask=REFERENCE_COMBINED_MASK_TSV,
                hard_json=REFERENCE_HARD_MASK_JSON,
                soft_json=REFERENCE_SOFT_MASK_JSON,
                dynamic_json=REFERENCE_DYNAMIC_MASK_JSON,
                combined_json=REFERENCE_COMBINED_MASK_JSON,
                summary_json=REFERENCE_MASK_SUMMARY_JSON
            log:
                project_path("logs", "reference_assets", "build_bin_mask.log")
            params:
                python_bin=config["biosoft"]["python"],
                script=SCRIPT_BUILD_BIN_ANNOTATIONS,
                script_action=SCRIPT_BUILD_BIN_ANNOTATIONS_ACTION,
                hard_n_fraction=HARD_MASK_N_FRACTION,
                soft_gc_low=SOFT_MASK_GC_LOW,
                soft_gc_high=SOFT_MASK_GC_HIGH,
                dynamic_z_frac=DYNAMIC_MASK_Z_FRAC,
                dynamic_median_abs_z=DYNAMIC_MASK_MEDIAN_ABS_Z
            threads: 1
            shell:
                r"""
                mkdir -p "$(dirname {output.hard_mask})" "$(dirname {log})"
                {params.python_bin:q} {params.script:q} {params.script_action:q} \
                    mask \
                    --annotation-tsvs {input.atomic_annotations:q} {input.analysis_annotations:q} {input.qc_annotations:q} \
                    --hard-mask-output {output.hard_mask:q} \
                    --soft-mask-output {output.soft_mask:q} \
                    --dynamic-mask-output {output.dynamic_mask:q} \
                    --combined-mask-output {output.combined_mask:q} \
                    --hard-mask-json-output {output.hard_json:q} \
                    --soft-mask-json-output {output.soft_json:q} \
                    --dynamic-mask-json-output {output.dynamic_json:q} \
                    --combined-mask-json-output {output.combined_json:q} \
                    --summary-json-output {output.summary_json:q} \
                    --hard-n-fraction {params.hard_n_fraction} \
                    --soft-gc-low {params.soft_gc_low} \
                    --soft-gc-high {params.soft_gc_high} \
                    --dynamic-z-frac-threshold {params.dynamic_z_frac} \
                    --dynamic-median-abs-z-threshold {params.dynamic_median_abs_z} \
                    >> {log:q} 2>&1
                """

if REFERENCE_PACKAGE_TARGET_FILES:
    rule build_frozen_reference_package:
        input:
            refs=REF_TARGET_FILES,
            atomic_bins=REFERENCE_ATOMIC_BINS,
            analysis_bins=REFERENCE_ANALYSIS_BINS,
            qc_bins=REFERENCE_QC_BINS,
            atomic_annotations=REFERENCE_ATOMIC_BIN_ANNOTATIONS,
            analysis_annotations=REFERENCE_ANALYSIS_BIN_ANNOTATIONS,
            qc_annotations=REFERENCE_QC_BIN_ANNOTATIONS,
            hard_mask=REFERENCE_HARD_MASK_TSV,
            soft_mask=REFERENCE_SOFT_MASK_TSV,
            dynamic_mask=REFERENCE_DYNAMIC_MASK_TSV,
            combined_mask=REFERENCE_COMBINED_MASK_TSV,
            metadata=RUN_METADATA
        output:
            manifest=REFERENCE_PACKAGE_MANIFEST,
            inventory=REFERENCE_PACKAGE_INVENTORY,
            readme=REFERENCE_PACKAGE_README,
            done=REFERENCE_PACKAGE_DONE
        log:
            project_path("logs", "reference_assets", "build_frozen_reference_package.log")
        params:
            python_bin=config["biosoft"]["python"],
            script=SCRIPT_BUILD_FROZEN_REFERENCE_PACKAGE,
            script_action=SCRIPT_BUILD_FROZEN_REFERENCE_PACKAGE_ACTION,
            package_name=REFERENCE_PACKAGE_NAME,
            package_dir=REFERENCE_PACKAGE_DIR
        threads: 1
        shell:
            r"""
            mkdir -p "{params.package_dir}" "$(dirname {log})"
            {params.python_bin:q} {params.script:q} {params.script_action:q} \
                --package-name {params.package_name:q} \
                --package-dir {params.package_dir:q} \
                --reference-files {input.refs:q} \
                --asset-files {input.atomic_bins:q} {input.analysis_bins:q} {input.qc_bins:q} {input.atomic_annotations:q} {input.analysis_annotations:q} {input.qc_annotations:q} {input.hard_mask:q} {input.soft_mask:q} {input.dynamic_mask:q} {input.combined_mask:q} \
                --metadata-files {input.metadata:q} \
                --manifest-output {output.manifest:q} \
                --inventory-output {output.inventory:q} \
                --readme-output {output.readme:q} \
                --done-output {output.done:q} \
                >> {log:q} 2>&1
            """
