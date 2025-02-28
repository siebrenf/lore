rule pigeon_sort:
    input:
        gff=rules.isoseq_collapse.output.gff,
    output:
        gff=expand("{isoseq_collapse_dir}/{{sample}}.gff", **config),
    log:
        expand("{isoseq_collapse_dir}/{{sample}}_sort.log", **config),
    threads: 1
    resources:
        mem_mb=500,
    shell:
        """
        pigeon sort \
            --log-file {log} \
            --log-level TRACE \
            -o {output.gff} \
            {input.gff} > {log} 2>&1
        """


rule pigeon_prepare:
    """
    Sort and index the genome annotation, (optional) CAGE peak, and (optional) intropolis files before classification. 
    This step ensures that all records for a given chromosome/scaffold are contiguous within the file. 
    Additionally, if a reference fasta is provided, the fai index will be generated.

    Sources:
      - https://isoseq.how/classification/workflow.html#prepare-reference-files-for-pigeon
    """
    input:
        ref_fasta=config["genome_fasta"],  # rules.genome.output.fasta,
        ref_gtf=config["gene_annotation_gtf"],  # rules.annotation.output.gtf,
    output:
        txt=expand("{pigeon_dir}/{genome}_prepared.txt", **config),
    log:
        expand("{pigeon_dir}/{genome}_prepare.log", **config),
    threads: 1
    resources:
        mem_mb=500,
    shell:
        """
        pigeon prepare \
            --log-file {log} \
            --log-level TRACE \
            {input.ref_fasta} \
            {input.ref_gtf}  > {log} 2>&1
            
        touch {output.txt}
        """


rule pigeon_classify:
    """
    Classify isoforms into categories.
     
    Sources:
      - https://isoseq.how/classification/workflow.html#classify-isoforms
    """
    input:
        gff=rules.pigeon_sort.output.gff,
        flcn=rules.isoseq_collapse.output.abundance,
        ref_gtf=config["gene_annotation_gtf"],  # rules.annotation.output.gtf,
        ref_fasta=config["genome_fasta"],  # rules.genome.output.fasta,
    output:
        gff=expand("{pigeon_dir}/{{sample}}.gff", **config),
    log:
        expand("{pigeon_dir}/{{sample}}_classify.log", **config),
    params:
        outdir=config["pigeon_dir"],
        # TODO: find the origins of polyA.txt & TSS.bed
        polya=f"--poly-a {config["poly-a"]}" if "poly-a" in config else "",
        cage_peak=f"--cage-peak {config["tss_bed"]}" if "tss_bed" in config else "",
    threads: 1
    resources:
        mem_mb=500,
    shell:
        """
        pigeon classify \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --out-dir {params.outdir} \
            {params.polya} \
            {params.cage_peak} \
            --flcn {input.flcn} \
            {input.gff} \
            {input.ref_gtf} \
            {input.ref_fasta} > {log} 2>&1
        """


rule pigeon_filter:
    """
    Filter isoforms from the classification output.

    Sources:
      - https://isoseq.how/classification/workflow.html#filter-isoforms
    """
    input:
        gff=rules.pigeon_sort.output.gff,
    output:
        txt=expand("{pigeon_dir}/{{sample}}.txt", **config),
    log:
        expand("{pigeon_dir}/{{sample}}_filter.log", **config),
    params:
        outdir=config["pigeon_dir"],
    threads: 1
    resources:
        mem_mb=500,
    shell:
        """
        pigeon filter \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            {input.gff} \
            {output.txt} > {log} 2>&1
        """
