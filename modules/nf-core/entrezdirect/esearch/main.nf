process ENTREZDIRECT_ESEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(term)
    val database

    output:
    tuple val(meta), path("*.xml") , emit: xml
    tuple val("${task.process}"), val('ENTREZDIRECT'), eval('esearch -version 2>&1'), emit: versions_esearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    esearch \\
        -db $database \\
        -query $term \\
        $args > ${prefix}.xml

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.xml

    """
}
