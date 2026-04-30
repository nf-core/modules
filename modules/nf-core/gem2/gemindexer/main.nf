process GEM2_GEMINDEXER {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gem2:20200110--h9ee0642_1':
        'quay.io/biocontainers/gem2:20200110--h9ee0642_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gem"), emit: index
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('gem2'), val("20200110"), emit: versions_gem2, topic: versions

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
        --mm-tmp-prefix ./tmp \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gem
    touch ${prefix}.log
    """
}
