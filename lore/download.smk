from snakemake.io import expand


rule qc_scripts:
    output:
        detail=expand("{scripts_dir}/report_detail.R", **config),
        summary=expand("{scripts_dir}/report_summary.R", **config),
        knee=expand("{scripts_dir}/plot_knees.py", **config),
    log:
        expand("{results_dir}/qc_scripts.log", **config),
    params:
        prefix="https://raw.githubusercontent.com/PacificBiosciences/barcoding/refs/heads/master/scripts/r",
    resources:
        parallel_downloads=1,
    shell:
        """
        wget \
            --output-file {log} \
            --output-document={output.detail} \
            {params.prefix}/report_detail.R
            
        wget \
            --output-file {log} \
            --output-document={output.summary} \
            {params.prefix}/report_summary.R

        wget \
            --output-file {log} \
            --output-document={output.knee} \
            https://downloads.pacbcloud.com/public/dataset/MAS-Seq/PLOT-scripts/plot_knees.py
        """


rule adapters:
    output:
        fasta=expand("{results_dir}/adapters.fasta", **config),
    log:
        expand("{results_dir}/adapters.log", **config),
    params:
        url=config["adapters"],
    resources:
        parallel_downloads=1,
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
    log:
        expand("{results_dir}/primers.log", **config),
    params:
        url=config["primers"],
    resources:
        parallel_downloads=1,
    shell:
        """
        wget \
            --output-file {log} \
            --output-document={output.fasta} \
            {params.url}
        """


rule barcodes:
    output:
        txt=expand("{results_dir}/barcodes.txt", **config),
    log:
        expand("{results_dir}/barcodes.log", **config),
    params:
        url=config["barcodes"],
    resources:
        parallel_downloads=1,
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


# rule genome:
#     output:
#         fasta=expand("{results_dir}/{{genome}}.fa",**config),
#     log: expand("{results_dir}/{{genome}}_fasta.log",**config),
#     params:
#         url=config["genome_fasta"],
#     resources:
#         parallel_downloads=1,
#     shell:
#         """
#         wget \
#             --output-file {log} \
#             --output-document={output.fasta} \
#             {params.url}
#         """
#
#
# rule annotation:
#     output:
#         gtf=expand("{results_dir}/{{genome}}.gtf",**config),
#     log: expand("{results_dir}/{{genome}}_annotation.log",**config),
#     params:
#         url=config["genome_annotation"],
#     resources:
#         parallel_downloads=1,
#     shell:
#         """
#         wget \
#             --output-file {log} \
#             --output-document={output.gtf} \
#             {params.url}
#         """
