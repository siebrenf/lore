from snakemake.io import expand, directory


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


rule multiqc:
    """
    Aggregate Quality Control data into a single MultiQC report.
    """
    input:
        files=set(
            [
                config["skera_dir"],
                config["lima_dir"],
                config["isoseq_refine_dir"],
                config["isoseq_correct_dir"],
                config["pbmm2_dir"],
                config["isoseq_collapse_dir"],
                config["pigeon_classify_dir"],
                config["pigeon_report_dir"],
                *[f'{config["qc_dir"]}/{sample}' for sample in config["samples"]],
            ]
        ),
    output:
        report=expand("{qc_dir}/multiqc_report.html", **config),
        data=directory(expand("{qc_dir}/multiqc_report_data", **config)),
    params:
        outdir=config["qc_dir"],
        schema=workflow.source_path("schemas/multiqc_config.yaml"),
    log:
        expand("{qc_dir}/multiqc.log", **config),
    benchmark:
        expand("{benchmark_dir}/multiqc.txt", **config)[0]
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc \
        {input.files} \
        --no-ai \
        --force \
        --module samtools \
        --module custom_content \
        --outdir {params.outdir} \
        --filename multiqc_report.html \
        --config {params.schema} \
        --verbose > {log} 2>&1
        """
