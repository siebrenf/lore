from snakemake.io import expand


rule bam2bigwig:
    """
    Convert a bam file into a bigwig file.
    """
    input:
        bam=rules.pbmm2_align.output.bam,
        bai=rules.pbmm2_align.output.bai,
    output:
        bw=expand("{pbmm2_dir}/{{sample}}.bw", **config),
    log:
        expand("{pbmm2_dir}/deeptools_bamcoverage_{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/deeptools_bamcoverage_{{sample}}.txt", **config)[0]
    conda:
        "envs/deeptools.yaml"
    threads: 16
    shell:
        """
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.bw} \
            --numberOfProcessors {threads} \
            --normalizeUsing BPM \
            --binSize 1 \
            --verbose \
            > {log} 2>&1
        """
