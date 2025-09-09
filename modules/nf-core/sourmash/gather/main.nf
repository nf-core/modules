process SOURMASH_GATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.14--hdfd78af_0':
        'biocontainers/sourmash:4.8.14--hdfd78af_0' }"

    input:
    tuple val(meta), path(signature)
    path(database)
    val save_unassigned
    val save_matches_sig
    val save_prefetch
    val save_prefetch_csv

    output:
    tuple val(meta), path('*.csv.gz')             , optional:true, emit: result
    tuple val(meta), path('*_unassigned.sig.zip') , optional:true, emit: unassigned
    tuple val(meta), path('*_matches.sig.zip')    , optional:true, emit: matches
    tuple val(meta), path('*_prefetch.sig.zip')   , optional:true, emit: prefetch
    tuple val(meta), path('*_prefetch.csv.gz')    , optional:true, emit: prefetchcsv
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unassigned  = save_unassigned   ? "--output-unassigned ${prefix}_unassigned.sig.zip" : ''
    def matches     = save_matches_sig  ? "--save-matches ${prefix}_matches.sig.zip"         : ''
    def prefetch    = save_prefetch     ? "--save-prefetch ${prefix}_prefetch.sig.zip"       : ''
    def prefetchcsv = save_prefetch_csv ? "--save-prefetch-csv ${prefix}_prefetch.csv.gz"    : ''
    """
    sourmash gather \\
        $args \\
        --output ${prefix}.csv.gz \\
        ${unassigned} \\
        ${matches} \\
        ${prefetch} \\
        ${prefetchcsv} \\
        ${signature} \\
        ${database}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.csv.gz
    touch ${prefix}_unassigned.sig.zip
    touch ${prefix}_matches.sig.zip
    touch ${prefix}_prefetch.sig.zip
    echo "" | gzip > ${prefix}_prefetch.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

}
