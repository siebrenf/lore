from snakemake.io import expand


rule pbmm2_index:
    """
    Generate an genomic index for pbmm2_align.
    """
    input:
        ref_fasta=config["genome_fasta"],  # rules.genome.output.fasta,
    output:
        mmi=expand("{pbmm2_dir}/{genome}.mmi", **config),
    log:
        expand("{pbmm2_dir}/{genome}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_mmi_{genome}.txt", **config)[0]
    threads: 3
    resources:
        mem_mb=12_000,
    shell:
        """
        pbmm2 index \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --preset ISOSEQ \
            {input.ref_fasta} \
            {output.mmi} > {log} 2>&1
        """


rule pbmm2_align:
    """
    Map reads to a reference genome
    
    Sources:
      - https://isoseq.how/classification/workflow.html#map-reads-to-a-reference-genome
    """
    input:
        bam=rules.isoseq_groupdedup.output.bam,
        mmi=rules.pbmm2_index.output.mmi,
    output:
        bam=expand("{pbmm2_dir}/{{sample}}.bam", **config),
    log:
        expand("{pbmm2_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_mma_{{sample}}.txt", **config)[0]
    threads: 12
    resources:
        mem_mb=20_000,
    shell:
        """
        pbmm2 align \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --preset ISOSEQ \
            --sort \
            --unmapped \
            {input.mmi} \
            {input.bam} \
            {output.bam} > {log} 2>&1
        """
