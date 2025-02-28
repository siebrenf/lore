from snakemake.io import expand


rule adapters:
    output:
        fasta=expand("{results_dir}/adapters.fasta", **config),
    log: expand("{results_dir}/adapters.log", **config),
    params:
        url=config["adapters"],
    shell:
        """
        wget \
            --output-file {log} \
            --output-document={output.fasta} \
            {params.url}
        """


rule primers:
    output:
        fasta=expand("{results_dir}/primers.fasta", **config),
    log: expand("{results_dir}/primers.log", **config),
    params:
        url=config["primers"],
    shell:
        """
        wget \
            --output-file {log} \
            --output-document={output.fasta} \
            {params.url}
        """


rule barcodes:
    output:
        txt=expand("{results_dir}/barcodes.txt",**config),
    log: expand("{results_dir}/barcodes.log",**config),
    params:
        url=config["barcodes"],
    shell:
        """
        wget \
            --output-file {log} \
            --output-document={output.txt} \
            {params.url}

        # unzip the barcodes, if needed
        if [[ {params.url} =~ \\.gz$ ]]; then
            mv {output.txt} {output.txt}.gz 
            gunzip {output.txt}.gz 
        fi
        """
