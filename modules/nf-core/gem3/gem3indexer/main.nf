process GEM3_GEM3INDEXER {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gem3-mapper:3.6.1--h9d449c0_12':
        'quay.io/biocontainers/gem3-mapper:3.6.1--h9d449c0_12' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gem") , emit: index
    tuple val(meta), path("*.info"), emit: info
    tuple val("${task.process}"), val('gem3-indexer'), val("3.6.1"), emit: versions_gem3indexer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gem-indexer \\
        -i ${fasta} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gem
    touch ${prefix}.info
    """
}
