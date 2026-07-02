process SUPPA_GENERATEEVENTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.4--pyhdfd78af_0':
        'quay.io/biocontainers/suppa:2.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gtf)
    val format
    val event_type
    val boundary_type
    val variability_threshold
    val pool_genes

    output:
    tuple val(meta), path("*.events*.${format}")    , emit: events
    tuple val(meta), path("*.events*.gtf")          , emit: events_gtf, optional: true
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def event_type_str = event_type instanceof List ? event_type.join(' ') : event_type
    def boundary_type_arg = boundary_type ?: 'S'
    def threshold = variability_threshold ?: 10
    def pg = pool_genes ? '-p' : ''
    """
    suppa.py \\
        generateEvents \\
        $args \\
        -f ${format} \\
        -e ${event_type_str} \\
        -b ${boundary_type_arg} \\
        -t ${threshold} \\
        ${pg} \\
        -i ${gtf} \\
        -o ${prefix}.events
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (format == "ioi") {
        """
        echo $args
        touch ${prefix}.events.${format}
        """
    } else {
        """
        echo $args
        touch ${prefix}.events.${format}
        touch ${prefix}.events.gtf
        """
    }
}
