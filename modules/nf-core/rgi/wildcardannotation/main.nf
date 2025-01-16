process RGI_WILDCARDANNOTATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_0':
        'biocontainers/rgi:6.0.3--pyha8f3691_0' }"

    input:
    path(wildcard)
    path(card_json)

    output:
    path("wildcard_database_processed") , emit: db
    env RGI_VERSION                     , emit: tool_version
    env DB_VERSION                      , emit: db_version
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    rgi wildcard_annotation \\
        -i $wildcard \\
        --card_json $card_json \\
        -v 
        $args

    mkdir wildcard_database_processed
    cp *.gz wildcard_database_processed
    cp Variants-Download-README.txt wildcard_database_processed

    RGI_VERSION=\$(rgi main --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
    END_VERSIONS
    """

    stub:
    """
    touch card.fasta
    touch card_all.fasta

    mkdir wildcard_database_processed
    cp *.gz wildcard_database_processed
    cp Variants-Download-README.txt wildcard_database_processed

    RGI_VERSION=\$(rgi main --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
    END_VERSIONS
    """
}
