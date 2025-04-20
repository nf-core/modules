process NAIL_SEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nail:0.3.0--h4349ce8_0':
        'biocontainers/nail:0.3.0--h4349ce8_0' }"

    input:
    tuple val(meta), path(query)
    path target
    val write_align

    output:
    tuple val(meta), path("${prefix}.txt"), emit: output
    tuple val(meta), path('results.tbl')  , emit: target_summary
    tuple val(meta), path("${prefix}.ali"), emit: alignments    , optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    alignment   = write_align ? "--ali-out ${prefix}.ali" : ''
    def VERSION = '0.3.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    nail search \\
        $args \\
        $alignment \\
        ${query} \\
        ${target} > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nail: $VERSION
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.txt
    ${write_align ? "touch ${prefix}.ali" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nail: $VERSION
    END_VERSIONS
    """
}
