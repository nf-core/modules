process RATTLE_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rattle:1.0--h5ca1c30_0' :
        'quay.io/biocontainers/rattle:1.0--h5ca1c30_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("clusters.out"), emit: clusters
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('rattle'), val("1.0"), emit: versions_rattle, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rattle \\
        cluster ${args} \\
        -t ${task.cpus} \\
        -i ${reads}
    """

    stub:
    """
    touch clusters.out
    """
}
