from snakemake.io import expand


rule rule_name:
    """
    docstring
    """
    input:
        key=rules.other_rule_name.output.key,
    output:
        key=expand("{config_key}/{{wildcard}}.bam", **config),
    log:
        "",
    params:
        flags="-a",
    threads: 1
    resources:
        mem_mb=40_000,
    shell:
        """
        touch \
        {params.flags} \
        {output.key} \
        > {log} 2>&1
        """
