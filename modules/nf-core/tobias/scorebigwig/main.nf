process TOBIAS_SCOREBIGWIG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tobias:0.17.5--py310h3479294_0':
        'quay.io/biocontainers/tobias:0.17.5--py310h3479294_0' }"

    input:
    tuple val(meta), path(signal), path(regions)

    output:
    tuple val(meta), path("*_footprints.bw"), emit: footprints
    tuple val("${task.process}"), val('tobias'), eval('TOBIAS --version'), topic: versions, emit: versions_tobias

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p .matplotlib
    export MPLCONFIGDIR="\${PWD}/.matplotlib"

    TOBIAS ScoreBigwig \\
        --signal $signal \\
        --regions $regions \\
        --output ${prefix}_footprints.bw \\
        --cores ${task.cpus} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_footprints.bw
    """
}
