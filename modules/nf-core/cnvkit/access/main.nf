process CNVKIT_ACCESS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/3935f8d7507f85fd20e073f0625bb6bdff8aa1b6044f4bbb720c9a8063ea0390/data'
:         'community.wave.seqera.io/library/cnvkit:0.9.14--288e98d6210b7304' }"

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
