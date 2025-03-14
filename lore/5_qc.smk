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
        txt=expand("{qc_dir}/lima_qc_detail_{{sample}}.txt", **config),
    params:
        outdir=config["qc_dir"],
    log:
        expand("{qc_dir}/lima_qc_detail_{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/lima_qc_detail_{{sample}}.txt", **config)[0]
    threads: 1
    resources:
        mem_mb=500,
    conda:
        "envs/lima_qc.yaml"
    shell:
        """
		CURDIR=$(pwd)
		
		cd {params.outdir}
        Rscript --vanilla {input.script} {input.report} > {log} 2>&1
		cd $CURDIR
		
		# TODO: rename output and use out official output
		touch {output.txt}
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
		expand("{qc_dir}/{{sample}}_summary_bad_adapter_ratio.png", **config),
		expand("{qc_dir}/{{sample}}_summary_hq_length_hist_2d.png", **config),
		expand("{qc_dir}/{{sample}}_summary_meanscore_vs_yield_hex.png", **config),
		expand("{qc_dir}/{{sample}}_summary_meanscore_vs_yield_hex_log10.png", **config),
		expand("{qc_dir}/{{sample}}_summary_meanscore_vs_yield_jitter.png", **config),
		expand("{qc_dir}/{{sample}}_summary_meanscore_vs_yield_jitter_log10.png", **config),
		expand("{qc_dir}/{{sample}}_summary_read_length_hist_2d.png", **config),
		expand("{qc_dir}/{{sample}}_summary_score_hist.png", **config),
		expand("{qc_dir}/{{sample}}_summary_score_hist_2d.png", **config),
		expand("{qc_dir}/{{sample}}_summary_yield_zmw.png", **config)
    params:
        outdir=config["qc_dir"],
    log:
        expand("{qc_dir}/lima_qc_summary_{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/lima_qc_summary_{{sample}}.txt", **config)[0]
    threads: 1
    resources:
        mem_mb=500,
    conda:
        "envs/lima_qc.yaml"
    shell:
        """
		CURDIR=$(pwd)
		
		cd {params.outdir}
        Rscript --vanilla {input.script} {input.report} > {log} 2>&1
		cd $CURDIR
		
		# rename the output to include the sample name
		declare -a List=(
			"summary_bad_adapter_ratio.png" 
			"summary_hq_length_hist_2d.png" 
			"summary_meanscore_vs_yield_hex.png"
			"summary_meanscore_vs_yield_hex_log10.png"
			"summary_meanscore_vs_yield_jitter.png"
			"summary_meanscore_vs_yield_jitter_log10.png"
			"summary_read_length_hist_2d.png"
			"summary_score_hist.png"
			"summary_score_hist_2d.png"
			"summary_yield_zmw.png"
        )
		for fname in "${List[@]}"; do
			mv {params.outdir}/$fname {params.outdir}/{wildcards.sample}_$fname
		done
        """
