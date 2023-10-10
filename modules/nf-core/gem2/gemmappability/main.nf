process GEM2_GEMMAPPABILITY {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::gem2=20200110"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gem2:20200110--h9ee0642_1':
        'biocontainers/gem2:20200110--h9ee0642_1' }"

    input:
    tuple val(meta), path(index)
    val(read_length)

    output:
    tuple val(meta), path("*.mappability")  , emit: map
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20200110' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gem-mappability \\
        -I ${index} \\
        -o ${prefix} \\
        -l ${read_length} \\
        -T ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem2: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20200110' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.mappability

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem2: $VERSION
    END_VERSIONS
    """
}
