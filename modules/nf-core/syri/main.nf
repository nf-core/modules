process SYRI {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/syri:1.7.0--py310hdbdd923_0':
        'biocontainers/syri:1.7.0--py310hdbdd923_0' }"

    input:
    tuple val(meta), path(infile)
    path(query_fasta)
    path(reference_fasta)
    val(file_type)

    output:
    tuple val(meta), path("*syri.out")      , emit: syri        , optional: true
    tuple val(meta), path("*.error.log")    , emit: error       , optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( ! ( file_type in [ 'T', 'S', 'B', 'P' ] ) ) { error "File type should be one of [ 'T', 'S', 'B', 'P' ]" }
    """
    syri \\
        -c $infile \\
        -q $query_fasta \\
        -r $reference_fasta \\
        -F $file_type \\
        $args \\
        --prefix $prefix \\
        2> >(tee "${prefix}.error.log" >&2) \\
        || echo "Errors from syri printed to ${prefix}.error.log"

    [ -f "${prefix}syri.out" ] \\
        && rm "${prefix}.error.log" \\
        || echo 'Syri failed and no syri.out file was created'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        syri: \$(syri --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}syri.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        syri: \$(syri --version)
    END_VERSIONS
    """
}
