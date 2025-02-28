from snakemake.io import expand


rule lima:
    """
    docstring
    """
    input:
        bam=expand("{ccs_dir}/{{sample}}.bam", **config),
        adapters=config["adapters"],  # ,
    output:expand("{lima_dir}/{{sample}}.fl.bam", **config),
    log: expand("{results_dir}/lima_{{sample}}.log", **config),
    threads: 8
    resources:
        mem_mb=7_000,
    shell:
        """
        lima --isoseq \
            --num-threads {threads} \
            --log-file {log} \
            {input.bam} {input.adapters} \
            {output}
        """


rule isoseq_tag:
    """
    docstring
    """
    input:
        bam=rules.lima.output
    output:
        bam=expand("{isoseq_tag_dir}/{{sample}}.tagged.bam", **config),
    log: expand("{results_dir}/isoseq_tag_{{sample}}.log", **config),
    params:
        design=config["design"],
		outdir=config["isoseq_tag_dir"],
    threads: 8
    resources:
        mem_mb=7_000,
    shell:
        """
        isoseq tag \
            --log-file {log} \
            --num-threads {threads} \
            --design {params.design} \
            {input.bam} \
            {output.bam}
        """
