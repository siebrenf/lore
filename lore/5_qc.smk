from snakemake.io import expand


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
        expand("{qc_dir}/{{sample}}/summary_meanscore_vs_yield_jitter_log10.png", **config),
        expand("{qc_dir}/{{sample}}/summary_read_length_hist_2d.png", **config),
        expand("{qc_dir}/{{sample}}/summary_score_hist.png", **config),
        expand("{qc_dir}/{{sample}}/summary_score_hist_2d.png", **config),
        expand("{qc_dir}/{{sample}}/summary_yield_zmw.png", **config)
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
