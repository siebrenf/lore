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
        xml=expand("{isoseq_tag_dir}/{{sample}}.transcriptset.xml",**config),
    log: expand("{isoseq_tag_dir}/{{sample}}.log", **config),
    params:
        design=config["design"],  # Barcoding design. Specifies which bases to use as cell/molecular barcodes
    threads: 8  # TODO: check
    resources:
        mem_mb=7_000,  # TODO: check
    shell:
        """
        isoseq tag \
            --log-file {log} \
            --num-threads {threads} \
            --design {params.design} \
            {input.bam} \
            {output.bam}
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
        adapters=rules.adapters.output.fasta,
    output:
        bam=expand("{isoseq_refine_dir}/{{sample}}.bam", **config),
        xml=expand("{isoseq_refine_dir}/{{sample}}.transcriptset.xml",**config),
    log: expand("{isoseq_refine_dir}/{{sample}}.log", **config),
    threads: 8  # TODO: check
    resources:
        mem_mb=7_000,  # TODO: check
    shell:
        """
        isoseq refine \
            --require-polya \
            --log-file {log} \
            --num-threads {threads} \
            {input.bam} \
            {input.adapters} \
            {output.bam}
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
    log: expand("{isoseq_correct_dir}/{{sample}}.log", **config),
    threads: 8  # TODO: check
    resources:
        mem_mb=7_000,  # TODO: check
    shell:
        """
        isoseq correct \
            --log-file {log} \
            --num-threads {threads} \
            --barcodes {input.barcodes} \
            {input.bam} \
            {output.bam}
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
        bam=expand("{dedup_dir}/{{sample}}.bam", **config),
        pbi=expand("{dedup_dir}/{{sample}}.bam.pbi", **config),
    log: expand("{dedup_dir}/{{sample}}.log", **config),
    threads: 8  # TODO: check
    resources:
        mem_mb=7_000,  # TODO: check
    shell:
        """
        isoseq groupdedup \
            --log-file {log} \
            --num-threads {threads} \
            {input.bam} \
            {output.bam}
        """
