from snakemake.io import expand, temp


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
        pbi=expand("{isoseq_tag_dir}/{{sample}}.bam.pbi", **config),
    log:
        expand("{isoseq_tag_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_tag_{{sample}}.txt", **config)[0]
    params:
        flags=config["isoseq_tag"],
    threads: 16
    resources:
        mem_mb=500,
    shell:
        """
        isoseq tag \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            {params.flags} \
            {input.bam} \
            {output.bam} \
            --verbose > {log} 2>&1
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
      - https://isoseq.how/umi/isoseq-correct.html
    """
    input:
        bam=rules.isoseq_tag.output.bam,
        primers=rules.primers.output.fasta,
    output:
        bam=expand("{isoseq_refine_dir}/{{sample}}.bam", **config),
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
    params:
        flags=config["isoseq_refine"],
    threads: 16
    resources:
        mem_mb=500,
    shell:
        """
        isoseq refine \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            {params.flags} \
            {input.bam} \
            {input.primers} \
            {output.bam} \
            --verbose > {log} 2>&1
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
        pbi=expand("{isoseq_correct_dir}/{{sample}}.bam.pbi", **config),
        json=expand("{isoseq_correct_dir}/{{sample}}.report.json", **config),
        pbi2=temp(
            expand("{isoseq_correct_dir}/{{sample}}_intermediate.bam.pbi", **config)
        ),
    log:
        expand("{isoseq_correct_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_cor_{{sample}}.txt", **config)[0]
    threads: 16
    resources:
        mem_mb=30_000,
    shell:
        """
        isoseq correct \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --barcodes {input.barcodes} \
            {input.bam} \
            {output.bam} \
            --verbose > {log} 2>&1
        """


rule isoseq_bcstats:
    """
    Generates stats for group barcodes and (optionally) molecular barcodes
    
    Sources:
      - https://isoseq.how/umi/isoseq-bcstats.html
    """
    input:
        bam=rules.isoseq_correct.output.bam,
    output:
        tsv=expand("{isoseq_correct_dir}/{{sample}}_bcstats.tsv", **config),
        json=expand("{isoseq_correct_dir}/{{sample}}_bcstats.json", **config),
    log:
        expand("{isoseq_correct_dir}/{{sample}}_bcstats.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_bcstats_{{sample}}.txt", **config)[0]
    params:
        flags=config["isoseq_bcstats"],
    threads: 16
    resources:
        mem_mb=20_000,
    shell:
        """
        isoseq bcstats \
            --log-file {log} \
            --log-level TRACE \
            --num-threads {threads} \
            --json {output.json} \
            -o {output.tsv} \
            {params.flags} \
            {input.bam} \
            --verbose > {log} 2>&1
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
        pbi=expand("{isoseq_dedup_dir}/{{sample}}.bam.pbi", **config),
        fasta=expand("{isoseq_dedup_dir}/{{sample}}.fasta", **config),  # TODO: use for separate analysis?
    log:
        expand("{isoseq_dedup_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/isoseq_groupdedup_{{sample}}.txt", **config)[0]
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
            {output.bam} \
            --verbose > {log} 2>&1
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
            {output.gff} \
            --verbose > {log} 2>&1
        """
