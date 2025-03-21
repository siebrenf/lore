from snakemake.io import expand, directory, unpack


rule lima_qc_detail:
    """
    Visualize various QC values of a lima.report file
     
    Sources:
      - https://lima.how/faq/qc.html
    """
    input:
        script=rules.qc_scripts.output.detail,
        report=rules.lima.output.report,
    output:
        # TODO: add other output here
        expand("{qc_dir}/{{sample}}/detail_yield_zmw.png", **config),
    params:
        outdir=config["qc_dir"],
    log:
        expand("{qc_dir}/lima_qc_detail_{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/lima_qc_detail_{{sample}}.txt", **config)[0]
    threads: 1
    resources:
        mem_mb=3_000,
    conda:
        "envs/lima_qc.yaml"
    shell:
        """
        CURDIR=$(pwd)
        mkdir -p {params.outdir}/{wildcards.sample}
        cd {params.outdir}/{wildcards.sample}
        Rscript --vanilla {input.script} {input.report} > {log} 2>&1
        cd $CURDIR
        """


rule lima_qc_summary:
    """
    Visualize high-plex data of a lima.report file
     
    Sources:
      - https://lima.how/faq/qc.html
    """
    input:
        script=rules.qc_scripts.output.summary,
        report=rules.lima.output.report,
    output:
        expand("{qc_dir}/{{sample}}/summary_bad_adapter_ratio.png", **config),
        expand("{qc_dir}/{{sample}}/summary_hq_length_hist_2d.png", **config),
        expand("{qc_dir}/{{sample}}/summary_meanscore_vs_yield_hex.png", **config),
        expand("{qc_dir}/{{sample}}/summary_meanscore_vs_yield_hex_log10.png", **config),
        expand("{qc_dir}/{{sample}}/summary_meanscore_vs_yield_jitter.png", **config),
        expand(
            "{qc_dir}/{{sample}}/summary_meanscore_vs_yield_jitter_log10.png", **config
        ),
        expand("{qc_dir}/{{sample}}/summary_read_length_hist_2d.png", **config),
        expand("{qc_dir}/{{sample}}/summary_score_hist.png", **config),
        expand("{qc_dir}/{{sample}}/summary_score_hist_2d.png", **config),
        expand("{qc_dir}/{{sample}}/summary_yield_zmw.png", **config),
    params:
        outdir=config["qc_dir"],
    log:
        expand("{qc_dir}/lima_qc_summary_{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/lima_qc_summary_{{sample}}.txt", **config)[0]
    threads: 1
    resources:
        mem_mb=8_000,
    conda:
        "envs/lima_qc.yaml"
    shell:
        """
        CURDIR=$(pwd)
        mkdir -p {params.outdir}/{wildcards.sample}
        cd {params.outdir}/{wildcards.sample}
        Rscript --vanilla {input.script} {input.report} > {log} 2>&1
        cd $CURDIR
        """


def reads_per_step_input(wildcards):
    files = {
        "skera": [],
        "lima": [],
        "refine": [],
        "correct": [],
    }
    for sample in config["samples"]:
        files["skera"].append(f'{config["skera_dir"]}/{sample}.summary.json')
        files["lima"].append(f'{config["lima_dir"]}/{sample}.lima.summary')
        files["refine"].append(
            f'{config["isoseq_refine_dir"]}/{sample}.filter_summary.report.json'
        )
        files["correct"].append(f'{config["isoseq_correct_dir"]}/{sample}.report.json')
    return files


rule reads_per_step_qc:
    """
    Generate a table with the number of each per processing step for each sample.
    """
    input:
        unpack(reads_per_step_input),
    output:
        expand("{qc_dir}/reads_per_step.tsv", **config),
    script:
        "scripts/reads_per_step.py"


rule isoseq_correct_qc:
    input:
        [
            f'{config["isoseq_correct_dir"]}/{sample}.report.json'
            for sample in config["samples"]
        ],
    output:
        expand("{qc_dir}/isoseq_correct.tsv", **config),
    script:
        "scripts/isoseq_correct.py"


rule samtools_stats:
    input:
        bam=expand("{{path}}/{{sample}}.bam", **config),
    output:
        tsv=expand("{{path}}/{{sample}}.stats.tsv", **config),
    log:
        expand("{{path}}/{{sample}}.stats.log", **config),
    threads: 4
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools stats --threads {threads} {input.bam} 1> {output.tsv} 2> {log}
        """


rule samtools_coverage:
    input:
        bam=expand("{pbmm2_dir}/{{sample}}.bam", **config),
    output:
        tsv=expand("{pbmm2_dir}/{{sample}}.coverage.tsv", **config),
    log:
        expand("{pbmm2_dir}/{{sample}}.coverage.log", **config),
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools coverage {input.bam} --output {output.tsv} 2> {log}
        """


rule mt_nuc_ratio_calculator:
    """
    Estimate the amount of nuclear and mitochondrial reads in a sample. 
    """
    input:
        bam=rules.pbmm2_align.output.bam,
    output:
        txt=expand("{pbmm2_dir}/{{sample}}.bam.mtnucratio", **config),
        json=expand("{pbmm2_dir}/{{sample}}.bam.mtnucratiomtnuc.json", **config),
    log:
        expand("{pbmm2_dir}/{{sample}}.bam.mtnucratio.log", **config),
    params:
        mitochondria=config["mitochondria"],
    conda:
        "envs/mtnucratio.yaml"
    shell:
        """
        mtnucratio {input.bam} {params.mitochondria} > {log} 2>&1
        """


rule multiqc_schema:
    output:
        schema=temp(expand("{qc_dir}/multiqc_config.yaml", **config)),
        # schema=expand("{qc_dir}/multiqc_config.yaml", **config),
    params:
        schema=workflow.source_path("schemas/multiqc_config.yaml"),
        config=config,
    log:
        expand("{qc_dir}/multiqc_config.log", **config),
    script:
        "scripts/multiqc_schema.py"


def qc_files(wildcards):
    files = {
        "samtools_stats": [],
        "samtools_coverage": [],
        "mtnucratio": [],
    }
    for sample in config["samples"]:
        # files["samtools_stats"].append(f'{config["skera_dir"]}/{sample}.stats.tsv')
        files["samtools_stats"].append(f'{config["pbmm2_dir"]}/{sample}.stats.tsv')
        files["samtools_coverage"].append(
            f'{config["pbmm2_dir"]}/{sample}.coverage.tsv'
        )
        files["mtnucratio"].append(
            f'{config["pbmm2_dir"]}/{sample}.bam.mtnucratiomtnuc.json'
        )
    return files


rule multiqc:
    """
    Aggregate Quality Control data into a single MultiQC report.
    """
    input:
        # aggregated files
        rules.reads_per_step_qc.output,
        rules.isoseq_correct_qc.output,
        # per-sample files
        unpack(qc_files),
        # everything in the QC directory
        dir=config["qc_dir"],
        # misc
        schema=rules.multiqc_schema.output.schema,
    output:
        report=expand("{qc_dir}/multiqc_report.html", **config),
        data=directory(expand("{qc_dir}/multiqc_report_data", **config)),
    params:
        dir=config["qc_dir"],
    log:
        expand("{qc_dir}/multiqc.log", **config),
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc \
        {input} \
        --outdir {params.dir} \
        --filename multiqc_report.html \
        --config {input.schema} \
        --no-ai \
        --no-version-check \
        --force \
        --dirs-depth 1 \
        --strict \
        --module samtools \
        --module mtnucratio \
        --module custom_content \
        --verbose > {log} 2>&1
        """


rule webshare:
    """
    Move files to a webshare directory.
    """
    input:
        report=rules.multiqc.output.report
    output:
        report=expand("{www_dir}/multiqc_report.html", **config),
    log:
        expand("{qc_dir}/webshare.log", **config),
    params:
        dir=config.get("www_dir", "."),
    shell:
        """
        rsync -ar {input.report} {output.report} > {log} 2>&1

        chmod -R 755 {params.dir} >> {log} 2>&1
        """
