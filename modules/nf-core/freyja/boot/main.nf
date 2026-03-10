process FREYJA_BOOT {
    tag "$meta.id"
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:2.0.3--pyhdfd78af_0' :
        'biocontainers/freyja:2.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(variants), path(depths)
    val repeats
    path barcodes
    path lineages_meta
    path lineages_topology

    output:
    tuple val(meta), path("*lineages.csv")  , emit: lineages
    tuple val(meta), path("*summarized.csv"), emit: summarized
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
        boot \\
        $args \\
        --nt $task.cpus \\
        --nb $repeats \\
        --output_base $prefix \\
        --barcodes $barcodes \\
        $meta_cmd \\
        $lineage_cmd \\
        $variants \\
        $depths

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_lineages.csv
    touch ${prefix}_summarized.csv

    """
}
