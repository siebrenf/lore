from snakemake.io import expand


rule lima:
    """
    Removal of primers and identification of barcodes is performed using lima. 
    
    Sources:
      - https://isoseq.how/umi/cli-workflow.html#step-2---primer-removal
    """
    input:
        bam=rules.skera.output.bam,
        adapters=rules.primers.output.fasta,
    output:
        bam=expand("{lima_dir}/{{sample}}.bam", **config),
        # unrequested files:
        pbi=expand("{lima_dir}/{{sample}}.bam.pbi", **config),
        xml=expand("{lima_dir}/{{sample}}.consensusreadset.xml", **config),
        json=expand("{lima_dir}/{{sample}}.json", **config),
        clips=expand("{lima_dir}/{{sample}}.lima.clips", **config),
        counts=expand("{lima_dir}/{{sample}}.lima.counts", **config),
        report=expand("{lima_dir}/{{sample}}.lima.report", **config),
        summary=expand("{lima_dir}/{{sample}}.lima.summary", **config),
    log:
        expand("{lima_dir}/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/lima_{{sample}}.txt", **config)[0]
    params:
        dir=config["lima_dir"],
    threads: 8
    resources:
        mem_mb=7_000,
    shell:
        """
        lima \
            --isoseq \
            --num-threads {threads} \
            --log-file {log} \
            --log-level TRACE \
            {input.bam} \
            {input.adapters} \
            {output.bam} > {log} 2>&1
            
        # remove the adapter direction from the file names
        for f1 in {params.dir}/{wildcards.sample}*; do
            f2="$(echo "$f1" | sed 's/5p--3p.//g; s/3p--5p.//g')"
            if [[ "$f1" != "$f2" ]]; then
                mv -f $f1 $f2
            fi
        done
        """
