from snakemake.io import expand


rule isoseq_tag:
    """
    Tags, such as UMIs and cell barcodes, have to be clipped from the reads and 
    associated with the reads for later deduplication.
    
    Sources:
      - https://isoseq.how/umi/cli-workflow.html#step-3---tag
    """
    input:
        bam=rules.lima.output.bam,
    output:
        bam=expand("{isoseq_tag_dir}/{{sample}}.bam", **config),
        # unrequested files:
        pbi=expand("{isoseq_tag_dir}/{{sample}}.bam.pbi", **config),
    log:
        expand("{isoseq_tag_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_tag_{{sample}}.txt", **config)[0]
    params:
        design=config["design"],  # Barcoding design. Specifies which bases to use as cell/molecular barcodes
    threads: 16
    resources:
        mem_mb=500,
    shell:
        """
        isoseq tag \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --design {params.design} \
            {input.bam} \
            {output.bam} > {log} 2>&1
        """


rule isoseq_refine:
    """
    Your data now contains full-length tagged reads, but still needs to be refined by:
      - Trimming of poly(A) tails
      - Unintended concatmer identification and removal 
        (note: if the library was constructed using the MAS-Seq method, 
        the reads should have already gone through skera 
        and is not expected to contain any more concatemers at this step)

    Sources:
      - https://isoseq.how/umi/cli-workflow.html#step-4---refine
    """
    input:
        bam=rules.isoseq_tag.output.bam,
        primers=rules.primers.output.fasta,
    output:
        bam=expand("{isoseq_refine_dir}/{{sample}}.bam", **config),
        # unrequested files:
        pbi=expand("{isoseq_refine_dir}/{{sample}}.bam.pbi", **config),
        xml=expand("{isoseq_refine_dir}/{{sample}}.consensusreadset.xml", **config),
        json=expand(
            "{isoseq_refine_dir}/{{sample}}.filter_summary.report.json", **config
        ),
        csv=expand("{isoseq_refine_dir}/{{sample}}.report.csv", **config),
    log:
        expand("{isoseq_refine_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_ref_{{sample}}.txt", **config)[0]
    threads: 16
    resources:
        mem_mb=500,
    shell:
        """
        isoseq refine \
            --require-polya \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            {input.bam} \
            {input.primers} \
            {output.bam} > {log} 2>&1
        """


rule isoseq_correct:
    """
    Identifies cell barcode errors and corrects them. 
    The tool uses a cell barcode whitelist to reassign erroneous barcodes based on edit distance. 
    Additionally, correct estimates which reads are likely to originate from a real cell and labels them using the rc tag.

    Common single-cell whitelists (e.g. 10x whitelist for 3’ kit) can be found in the MAS-Seq dataset. 
    These are the reverse complement of the 10x single-cell whitelists.

    Sources:
      - https://isoseq.how/umi/cli-workflow.html#step-5---cell-barcode-correction-and-real-cell-identification
      - https://isoseq.how/umi/isoseq-correct.html#barcode-correction-documentation
    """
    input:
        bam=rules.isoseq_refine.output.bam,
        barcodes=rules.barcodes.output.txt,
    output:
        bam=expand("{isoseq_correct_dir}/{{sample}}.bam", **config),
        # unrequested files:
        pbi=expand("{isoseq_correct_dir}/{{sample}}.bam.pbi", **config),
        json=expand("{isoseq_correct_dir}/{{sample}}.report.json", **config),
        pbi2=expand("{isoseq_correct_dir}/{{sample}}_intermediate.bam.pbi", **config),
    log:
        expand("{isoseq_correct_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_cor_{{sample}}.txt", **config)[0]
    threads: 16  # TODO: check
    resources:
        mem_mb=40_000,  # TODO: check
    shell:
        """
        isoseq correct \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --barcodes {input.barcodes} \
            {input.bam} \
            {output.bam} > {log} 2>&1
        """


rule isoseq_groupdedup:
    """
    Performs PCR deduplication via clustering by UMI and cell barcodes (if available).
    Generates one consensus sequence per founder molecule, using a QV guided consensus.

    Sources:
      - https://isoseq.how/umi/cli-workflow.html#step-6---deduplication
      - https://isoseq.how/umi/dedup-faq.html
    """
    input:
        bam=rules.isoseq_correct.output.bam,
    output:
        bam=expand("{isoseq_dedup_dir}/{{sample}}.bam", **config),
        # unrequested files:
        pbi=expand("{isoseq_dedup_dir}/{{sample}}.bam.pbi", **config),
        fasta=expand("{isoseq_dedup_dir}/{{sample}}.fasta", **config),
    log:
        expand("{isoseq_dedup_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_gdd_{{sample}}.txt", **config)[0]
    threads: 8
    resources:
        mem_mb=8_000,
    shell:
        """
        isoseq groupdedup \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            {input.bam} \
            {output.bam} > {log} 2>&1
        """


rule isoseq_collapse:
    """
    Collapse redundant transcripts into unique isoforms based on exonic structures

    Sources:
      - https://isoseq.how/classification/workflow.html#collapse-into-unique-isoforms
    """
    input:
        bam=expand("{pbmm2_dir}/{{sample}}.bam", **config),  # bam=rules.pbmm2_align.output.bam,
    output:
        gff=expand("{isoseq_collapse_dir}/{{sample}}.unsorted.gff", **config),
        abundance=expand(
            "{isoseq_collapse_dir}/{{sample}}.unsorted.abundance.txt", **config
        ),
        # unrequested files:
        fasta=expand("{isoseq_collapse_dir}/{{sample}}.unsorted.fasta", **config),
        fastq=expand("{isoseq_collapse_dir}/{{sample}}.unsorted.fastq", **config),
        flnc_count=expand(
            "{isoseq_collapse_dir}/{{sample}}.unsorted.flnc_count.txt", **config
        ),
        group=expand("{isoseq_collapse_dir}/{{sample}}.unsorted.group.txt", **config),
        stats=expand(
            "{isoseq_collapse_dir}/{{sample}}.unsorted.read_stat.txt", **config
        ),
        json=expand("{isoseq_collapse_dir}/{{sample}}.unsorted.report.json", **config),
    log:
        expand("{isoseq_collapse_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_col_{{sample}}.txt", **config)[0]
    threads: 8
    resources:
        mem_mb=14_000,
    shell:
        """
        isoseq collapse \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            {input.bam} \
            {output.gff} > {log} 2>&1
        """
