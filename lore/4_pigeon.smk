rule pigeon_sort:
    input:
        gff=rules.isoseq_collapse.output.gff,
    output:
        gff=expand("{isoseq_collapse_dir}/{{sample}}.gff", **config),
    log:
        expand("{isoseq_collapse_dir}/{{sample}}_sort.log", **config),
    benchmark:
        expand("{benchmark_dir}/pigeon_srt_{{sample}}.txt", **config)[0]
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
        ref_gtf=config["gene_annotation_gtf"].replace(".gtf", ".sorted.gtf"),
        ref_pgi=config["gene_annotation_gtf"].replace(".gtf", ".sorted.gtf.pgi"),
    log:
        expand("{pigeon_dir}/{genome}_prepare.log", **config),
    benchmark:
        expand("{benchmark_dir}/pigeon_pre_{genome}.txt", **config)[0]
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
        ref_gtf=rules.pigeon_prepare.output.ref_gtf,
        ref_pgi=rules.pigeon_prepare.output.ref_pgi,
        ref_fasta=config["genome_fasta"],
    output:
        json=expand("{pigeon_dir}/{{sample}}.report.json", **config),
        summary=expand("{pigeon_dir}/{{sample}}.summary.txt", **config),
        classification=expand("{pigeon_dir}/{{sample}}_classification.txt", **config),
        junctions=expand("{pigeon_dir}/{{sample}}_junctions.txt", **config),
    log:
        expand("{pigeon_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/pigeon_cls_{{sample}}.txt", **config)[0]
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
            --flnc {input.flcn} \
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
        classification=rules.pigeon_classify.output.classification,
        gff=rules.pigeon_sort.output.gff,
    output:
        json=expand(
            "{pigeon_dir}/{{sample}}_classification.filtered.report.json", **config
        ),
        summary=expand(
            "{pigeon_dir}/{{sample}}_classification.filtered.summary.txt", **config
        ),
        classification=expand(
            "{pigeon_dir}/{{sample}}_classification.filtered_lite_classification.txt",
            **config,
        ),
        junctions=expand(
            "{pigeon_dir}/{{sample}}_classification.filtered_lite_junctions.txt",
            **config,
        ),
        reasons=expand(
            "{pigeon_dir}/{{sample}}_classification.filtered_lite_reasons.txt",
            **config,
        ),
    log:
        expand("{pigeon_dir}/{{sample}}_filter.log", **config),
    benchmark:
        expand("{benchmark_dir}/pigeon_flt_{{sample}}.txt", **config)[0]
    threads: 1
    resources:
        mem_mb=500,
    shell:
        """
        pigeon filter \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --isoforms {input.gff} \
            {input.classification} > {log} 2>&1
        """


rule pigeon_make_seurat:
    """
    Generate files for Seurat/Scanpy

    Sources:
      - https://isoseq.how/classification/workflow.html#report-gene-saturation
    """
    input:
        classification=rules.pigeon_filter.output.classification,
        fasta=rules.isoseq_groupdedup.output.fasta,
        group=rules.isoseq_collapse.output.group,
    output:
        genes_matrix=expand("{seurat_dir}/{{sample}}/genes_seurat/matrix.mtx", **config),
        isoforms_matrix=expand(
            "{seurat_dir}/{{sample}}/isoforms_seurat/matrix.mtx", **config
        ),
        info=expand("{seurat_dir}/{{sample}}/{{sample}}.info.csv", **config),
    log:
        expand("{seurat_dir}/{{sample}}_make_seurat.log", **config),
    benchmark:
        expand("{benchmark_dir}/pigeon_mks_{{sample}}.txt", **config)[0]
    params:
        outdir=expand("{seurat_dir}/{{sample}}", **config),
    threads: 4  # TODO: check
    resources:
        mem_mb=2_000,  # TODO: check
    shell:
        """
        pigeon make-seurat \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --out-dir {params.outdir} \
            --dedup {input.fasta} \
            --group {input.group} \
            --keep-novel-genes \
            --out-prefix {wildcards.sample} \
            {input.classification} > {log} 2>&1
        """


rule pigeon_report:
    """
    Gene and isoform- level saturation can be determined by subsampling the 
    classification output and determining the number of unique 
    genes / isoforms at each subsample size.

    Sources:
      - https://isoseq.how/classification/workflow.html#report-gene-saturation
    """
    input:
        classification=rules.pigeon_filter.output.classification,
    output:
        report=expand("{pigeon_dir}/{{sample}}_saturation.txt", **config),
    log:
        expand("{pigeon_dir}/{{sample}}_report.log", **config),
    benchmark:
        expand("{benchmark_dir}/pigeon_rep_{{sample}}.txt", **config)[0]
    threads: 2  # TODO: check
    resources:
        mem_mb=2_000,  # TODO: check
    shell:
        """
        pigeon report \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            {input.classification} \
            {output.report} > {log} 2>&1
        """
