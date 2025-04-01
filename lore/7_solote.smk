from snakemake.io import expand


rule dfam:
    """
    Downloads a TE database partition required to run RepeatMasker.
    Check out the sources to determine the required partition.
    
    Sources:
      - https://www.dfam.org/releases/current/families/FamDB/README.txt
      - https://github.com/Dfam-consortium/RepeatMasker
    """
    output:
        md5=expand("{solote_dir}/{dfam_partition}.h5.gz.md5",**config),
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
        echo "Source ${{URL}}" > {log}
        
        DB=${{CONDA_PREFIX}}/share/RepeatMasker/Libraries/famdb/{params.partition}.h5.gz
        echo "Targer ${{DB}}" >> {log}
        
        # wget \
        #     --quiet \
        #     --output-document=${{DB}} \
        #     ${{URL}}
        
        wget \
            --quiet \
            --output-document={output.md5} \
            ${{URL}}.md5
        MD5=$(md5sum ${{DB}} | awk '{{print $1}}')
        SUCCESS=$(grep --count -F $MD5 {output.md5})
        
        # Add md5sum check to the log file
        echo "MD5sums:" >> {log}
        md5sum ${{DB}} | awk '{{print $1}}' >> {log}
        cat {output.md5} | awk '{{print $1}}' >> {log}
        echo "" >> {log}
        
        # Use the md5sum check to validate the download
        if [ ! $SUCCESS == 1 ]; then
            exit 1
        fi
        
        gunzip ${{DB}}
        """


rule repeatmasker:
    """
    Map repeats in the genome
    """
    input:
        genome=config["genome_fasta"],
        db=rules.dfam.output,
    output:
        fasta=expand("{solote_dir}/{genome}.out.fa",**config),
    log:
        expand("{solote_dir}/{genome}_repeatmasker.log",**config),
    benchmark:
        expand("{benchmark_dir}/repeatmasker_{genome}.txt", **config)[0]
    params:
        outdir=config["solote_dir"],
        species=config["species"],
        flags=config["repeatmasker"],
    threads: 4
    conda:
        "envs/repeatmasker.yaml"
    shell:
        """
        RepeatMasker \
          -species "{params.species}" \
          -dir {params.outdir} \
          -pa {threads} \
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
#
#     List all precompiled databases with:
#       python ${CONDA_PREFIX}/bin/SoloTE_RepeatMasker_to_BED.py --list
#
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
#
#         python \
#           ${{CONDA_PREFIX}}/bin/SoloTE_RepeatMasker_to_BED.py \
#           --genome {params.genome} > {log} 2>&1
#
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
        # bed=rules.solote_database.output.bed,
        bed=rules.solote_custom_database.output.bed,
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
