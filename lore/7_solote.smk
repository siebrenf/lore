from snakemake.io import expand


rule solote_database:
    """
    Download a precompiled RepeatMasker .fa.out.gz file.  
    RepeatMasker was run with the -s (sensitive) setting.
    SoloTE converts the .fa.out.gz to a BED file.

    List all precompiled databases with:
      python ${CONDA_PREFIX}/bin/SoloTE_RepeatMasker_to_BED.py --list

    Sources:
      - https://github.com/bvaldebenitom/SoloTE
      - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
    """
    output:
        bed=expand("{solote_dir}/{solote_genome}_rmsk.bed", **config),
    log:
        expand("{solote_dir}/{solote_genome}_solote_database.log",**config),
    params:
        outdir=config["solote_dir"],
        genome=config["solote_genome"],
    conda:
        "envs/solote.yaml"
    shell:
        """
        CURDIR=$(pwd)
        cd {params.outdir}
        
        python \
          ${{CONDA_PREFIX}}/bin/SoloTE_RepeatMasker_to_BED.py \
          --genome {params.genome} > {log} 2>&1
          
        cd $CURDIR
        """


# python ${CONDA_PREFIX}/bin/SoloTE_pipeline.py --help
rule solote:
    """
    Analysis of transposable elements in single-cell RNA-Seq data using locus-specific expression

    Sources:
      - https://github.com/bvaldebenitom/SoloTE
    """
    input:
        bam=rules.pbmm2_align.output.bam,
        bed=rules.solote_database.output.bed,
    output:
        bed=expand("{solote_dir}/{{sample}}_does_not_exist.bam",**config),
    log:
        expand("{solote_dir}/{{sample}}.log",**config),
    benchmark:
        expand("{benchmark_dir}/solote_{{sample}}.txt", **config)[0]
    params:
        outdir=config["solote_dir"],
    threads: 60
    conda:
        "envs/solote.yaml"
    shell:
        """    
        python \
          ${{CONDA_PREFIX}}/bin/SoloTE_pipeline.py \
          --bam {input.bam} \
          --teannotation {input.bed} \
          --outputdir {params.outdir} \
          --outputprefix {wildcards.sample} \
          --threads {threads} \
          > {log} 2>&1
        """