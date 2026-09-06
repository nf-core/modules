process GEM2_GEM2BEDMAPPABILITY {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gem2:20200110--h9ee0642_1':
        'quay.io/biocontainers/gem2:20200110--h9ee0642_1' }"

    input:
    tuple val(meta) , path(map)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bg")   , emit: bedgraph
    tuple val(meta), path("*.sizes"), emit: sizes
    tuple val("${task.process}"), val('gem2'), val("20200110"), emit: versions_gem2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    gem-2-bed mappability \\
        --input ${map} \\
        --index ${index} \\
        --output ${prefix}
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bg
    touch ${prefix}.sizes
    """
}
