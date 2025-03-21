from snakemake.io import expand, temp


rule skera:
    """
    Deconcatenate Kinnex PacBio HiFi reads to produce S-reads that represent the 
    original cDNA molecules. 
    
    Sources: 
      - https://skera.how
      - https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/
    """
    input:
        bam=expand("{ccs_dir}/{{sample}}.bam", **config),
        adapters=rules.adapters.output.fasta,
    output:
        bam=expand("{skera_dir}/{{sample}}.bam", **config),
        pbi=expand("{skera_dir}/{{sample}}.bam.pbi", **config),
        lig=expand("{skera_dir}/{{sample}}.ligations.csv", **config),
        bam_np=temp(expand("{skera_dir}/{{sample}}.non_passing.bam", **config)),
        pbi_np=temp(expand("{skera_dir}/{{sample}}.non_passing.bam.pbi", **config)),
        rl=expand("{skera_dir}/{{sample}}.read_lengths.csv", **config),
        csv=expand("{skera_dir}/{{sample}}.summary.csv", **config),
        json=expand("{skera_dir}/{{sample}}.summary.json", **config),
    log:
        expand("{skera_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/skera_{{sample}}.txt", **config)[0]
    threads: 48
    resources:
        mem_mb=2_000,
    conda:
        "envs/skera.yaml"
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
