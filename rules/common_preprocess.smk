rule collect_run_metadata:
    output:
        RUN_METADATA
    log:
        project_path("logs", "metadata", "collect_run_metadata.log")
    params:
        python_bin=config["biosoft"]["python"],
        fastp=config["biosoft"]["fastp"],
        bwa=config["biosoft"]["bwa"],
        samtools=config["biosoft"]["samtools"],
        wise=config["biosoft"]["WisecondorX"],
        project_root=str(Path.cwd())
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        {params.python_bin:q} {SCRIPT_COLLECT_RUN_METADATA:q} \
            --output {output:q} \
            --project-root {params.project_root:q} \
            --fastp {params.fastp:q} \
            --bwa {params.bwa:q} \
            --samtools {params.samtools:q} \
            --wisecondorx {params.wise:q} \
            --python-bin {params.python_bin:q} \
            --log {log:q}
        """


rule fastp_bwa:
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2=lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        clean_r1=FASTP_R1,
        clean_r2=FASTP_R2,
        html=FASTP_HTML,
        json=FASTP_JSON,
        bam=SORTED_BAM,
        bai=SORTED_BAI
    log:
        fastp=project_path("logs", "fastp", "{sample}.log"),
        bwa=project_path("logs", "bwa", "{sample}.log")
    threads: 8
    params:
        fastp=config["biosoft"]["fastp"],
        bwa=config["biosoft"]["bwa"],
        samtools=config["biosoft"]["samtools"],
        reference=config["core"]["reference_genome"]
    shell:
        r"""
        mkdir -p "$(dirname {output.clean_r1})" "$(dirname {output.bam})" "$(dirname {log.fastp})" "$(dirname {log.bwa})"
        (
            echo "=== PIPELINE AUDIT ==="
            echo "git_branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || true)"
            echo "git_commit: $(git rev-parse HEAD 2>/dev/null || true)"
            echo "fastp_version: $({params.fastp:q} --version 2>&1 | head -n 1 || true)"
            echo "=== COMMAND ==="
        ) > {log.fastp:q}
        {params.fastp:q} \
            -i {input.r1:q} -I {input.r2:q} \
            -o {output.clean_r1:q} -O {output.clean_r2:q} \
            -h {output.html:q} -j {output.json:q} \
            -w {threads} \
            >> {log.fastp:q} 2>&1
        (
            echo "=== PIPELINE AUDIT ==="
            echo "git_branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || true)"
            echo "git_commit: $(git rev-parse HEAD 2>/dev/null || true)"
            echo "bwa_version: $({params.bwa:q} 2>&1 | head -n 1 || true)"
            echo "samtools_version: $({params.samtools:q} --version 2>&1 | head -n 1 || true)"
            echo "=== COMMAND ==="
        ) > {log.bwa:q}
        {params.bwa:q} mem -t {threads} {params.reference:q} {output.clean_r1:q} {output.clean_r2:q} 2>> {log.bwa:q} \
            | {params.samtools:q} sort -@ {threads} -o {output.bam:q}
        {params.samtools:q} index {output.bam:q}
        """
