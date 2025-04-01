from snakemake.io import expand


# ruleorder: solote_precompiled_database > solote_custom_database


rule dfam:
    """
    Downloads a TE database partition required to run RepeatMasker.
    Check out the sources to determine the required partition.
    
    Sources:
      - https://www.dfam.org/releases/current/families/FamDB/README.txt
      - https://github.com/Dfam-consortium/RepeatMasker
    """
    output:
        proof=expand("{solote_dir}/{dfam_partition}.done",**config)
    log:
        expand("{solote_dir}/{dfam_partition}.log",**config),
    params:
        outdir=config["solote_dir"],
        partition=config["dfam_partition"],
    conda:
        "envs/repeatmasker.yaml"
    shell:
        """
        URL=https://www.dfam.org/releases/current/families/FamDB/{params.partition}.h5.gz
        echo "Source: ${{URL}}" > {log}
        
        DB=${{CONDA_PREFIX}}/share/RepeatMasker/Libraries/famdb/{params.partition}.h5
        echo "Target: ${{DB}}" >> {log}

        if [ ! -f ${{DB}} ]; then
            wget \
                --quiet \
                --output-document=${{DB}}.gz \
                ${{URL}}
            
            gunzip ${{DB}}.gz
        fi

        touch {output.proof}
        """


rule repeatmasker:
    """
    Map repeats in the genome
    """
    input:
        genome=config["genome_fasta"],
        proof=rules.dfam.output,
    output:
        # TODO: add additional output
        fasta=expand("{solote_dir}/{genome}.out.fa",**config),
    log:
        expand("{solote_dir}/{genome}_repeatmasker.log",**config),
    benchmark:
        expand("{benchmark_dir}/repeatmasker_{genome}.txt", **config)[0]
    params:
        outdir=config["solote_dir"],
        species=config["species"],
        flags=config["repeatmasker"],
    threads: 124
    conda:
        "envs/repeatmasker.yaml"
    shell:
        """
        # Each parallel job uses 4 threads
        let PA={threads}/4

        RepeatMasker \
          -species "{params.species}" \
          -dir {params.outdir} \
          -pa ${{PA}} \
          {params.flags} \
          {input.genome} \
          > {log} 2>&1
        """


rule solote_custom_database:
    """
    Convert RepeatMasker output to a SoloTE input BED file.
    """
    input:
        rules.repeatmasker.output.fasta
    output:
        bed=expand("{solote_dir}/{genome}_rmsk.bed", **config),
    script:
        "scripts/repeatmasker2bed.py"


# rule solote_precompiled_database:
#     """
#     Download a precompiled RepeatMasker .fa.out.gz file.
#     RepeatMasker was run with the -s (sensitive) setting.
#     SoloTE converts the .fa.out.gz to a BED file.

#     List all precompiled databases with:
#       python ${CONDA_PREFIX}/bin/SoloTE_RepeatMasker_to_BED.py --list

#     Sources:
#       - https://github.com/bvaldebenitom/SoloTE
#       - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
#     """
#     output:
#         bed=expand("{solote_dir}/{solote_genome}_rmsk.bed", **config),
#     log:
#         expand("{solote_dir}/{solote_genome}_solote_database.log",**config),
#     params:
#         outdir=config["solote_dir"],
#         genome=config["solote_genome"],
#     conda:
#         "envs/solote.yaml"
#     shell:
#         """
#         CURDIR=$(pwd)
#         cd {params.outdir}

#         python \
#           ${{CONDA_PREFIX}}/bin/SoloTE_RepeatMasker_to_BED.py \
#           --genome {params.genome} > {log} 2>&1

#         cd $CURDIR
#         """


rule solote:
    """
    Analysis of transposable elements in single-cell RNA-Seq data using locus-specific expression

    Sources:
      - https://github.com/bvaldebenitom/SoloTE
      - python ${CONDA_PREFIX}/bin/SoloTE_pipeline.py --help
    """
    input:
        bam=rules.pbmm2_align.output.bam,
        bed=rules.solote_custom_database.output.bed,
    output:
        bed=expand("{solote_dir}/{{sample}}_does_not_exist.bam",**config),  # TODO: replace with actual output
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
