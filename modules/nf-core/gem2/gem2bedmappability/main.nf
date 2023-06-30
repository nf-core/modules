process GEM2_GEM2BEDMAPPABILITY {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::gem2=20200110"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gem2:20200110--h9ee0642_1':
        'biocontainers/gem2:20200110--h9ee0642_1' }"

    input:
    tuple val(meta) , path(map)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bg")   , emit: bedgraph
    tuple val(meta), path("*.sizes"), emit: sizes
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20200110' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gem-2-bed mappability \\
        --input ${map} \\
        --index ${index} \\
        --output ${prefix}

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
    touch ${prefix}.bg
    touch ${prefix}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem2: $VERSION
    END_VERSIONS
    """
}
