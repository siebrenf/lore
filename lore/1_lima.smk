from snakemake.io import expand


rule lima:
    """
    Removal of primers and identification of barcodes is performed using lima. 
    
    Sources:
      - https://isoseq.how/umi/cli-workflow.html#step-2---primer-removal
    """
    input:
        bam=rules.skera.output.bam,
        adapters=rules.primers.output.fasta,
    output:
        bam=expand("{lima_dir}/{{sample}}.bam", **config),
        build=expand("{lima_dir}/{{sample}}.bam.pbi.build", **config),
        clips=expand("{lima_dir}/{{sample}}.lima.clips", **config),
        report=expand("{lima_dir}/{{sample}}.lima.report", **config),
    log: expand("{lima_dir}/{{sample}}.log", **config),
    threads: 8
    resources:
        mem_mb=7_000,
    shell:
        """
        lima 
            --isoseq \
            --num-threads {threads} \
            --log-file {log} \
            {input.bam} \
            {input.adapters} \
            {output.bam}
        """
