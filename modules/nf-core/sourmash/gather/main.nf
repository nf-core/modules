process SOURMASH_GATHER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::sourmash=4.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.5.0--hdfd78af_0':
        'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(signature)
    path(database)

    output:
    tuple val(meta), path('*.csv')           , optional:true, emit: result
    tuple val(meta), path('*_unassigned.sig'), optional:true, emit: unassigned
    tuple val(meta), path('*_matches.sig')   , optional:true, emit: matches
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--ksize 31'
    def prefix = task.ext.prefix ?: "${meta.id}"
    result     = "${prefix}.csv"
    matches    = "${prefix}_matches.sig"
    unassigned = "${prefix}_unassigned.sig"

    """
    sourmash gather \\
        $args \\
        --output '${result}' \\
        --output-unassigned '${unassigned}' \\
        --save-matches '${matches}' \\
        '${signature}' \\
        '${database}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
