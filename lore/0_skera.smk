from snakemake.io import expand


rule skera:
    """
    Deconcatenate HiFi reads to produce S-reads that represent the original cDNA molecules. 
    
    Sources: 
      - https://skera.how
      - https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/
    """
    input:
        bam=expand("{ccs_dir}/{{sample}}.bam", **config),
        adapters=rules.adapters.output.fasta,
    output:
        bam=expand("{sreads_dir}/{{sample}}.bam", **config),
        # unrequested files:
        pbi=expand("{sreads_dir}/{{sample}}.bam.pbi", **config),
        lig=expand("{sreads_dir}/{{sample}}.ligations.csv", **config),
        bamnp=expand("{sreads_dir}/{{sample}}.non_passing.bam", **config),
        pbinp=expand("{sreads_dir}/{{sample}}.non_passing.bam.pbi", **config),
        rl=expand("{sreads_dir}/{{sample}}.read_lengths.csv", **config),
        csv=expand("{sreads_dir}/{{sample}}.summary.csv", **config),
        json=expand("{sreads_dir}/{{sample}}.summary.json", **config),
    log:
        expand("{sreads_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/skera_{{sample}}.txt", **config)[0]
    threads: 24  # TODO: check
    resources:
        mem_mb=7_000,  # TODO: check
    shell:
        """
        skera split \
            --num-threads {threads} \
            --log-file {log} \
            --log-level TRACE \
            {input.bam} \
            {input.adapters} \
            {output.bam} > {log} 2>&1
        """
