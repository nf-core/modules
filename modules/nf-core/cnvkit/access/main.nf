process CNVKIT_ACCESS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cnvkit:0.9.12--pyhdfd78af_0'
        : 'biocontainers/cnvkit:0.9.12--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(exclude_bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), topic: versions, emit: versions_cnvkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_cmd = exclude_bed.collect { bed_file -> "-x ${bed_file}" }.join(" ")
    """
    cnvkit.py \\
        access \\
        ${fasta} \\
        ${exclude_cmd} \\
        ${args} \\
        --output ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    """
}
