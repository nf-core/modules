process FREYJA_DEMIX {
    tag "$meta.id"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:2.0.3--pyhdfd78af_0' :
        'biocontainers/freyja:2.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(variants), path(depths)
    path barcodes
    path lineages_meta
    path lineages_topology

    output:
    tuple val(meta), path("*.tsv"), emit: demix
    tuple val("${task.process}"), val('freyja'), eval("freyja --version | sed 's/.* //'"), topic: versions, emit: versions_freyja

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def meta_cmd = lineages_meta ? "--meta $lineages_meta" : ''
    def lineage_cmd = lineages_topology ? "--lineageyml $lineages_topology" : ''
    """
    freyja \\
        demix \\
        $args \\
        --output ${prefix}.tsv \\
        --barcodes $barcodes \\
        $lineage_cmd \\
        $meta_cmd \\
        $variants \\
        $depths
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    """

}
