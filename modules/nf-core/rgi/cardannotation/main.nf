process RGI_CARDANNOTATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_0':
        'biocontainers/rgi:6.0.3--pyha8f3691_0' }"

    input:
    path(card)

    output:
    path("output_dir")  , emit: db
    env RGI_VERSION     , emit: tool_version
    env DB_VERSION      , emit: db_version
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    rgi card_annotation -i $card/card.json \\
    $args

    mkdir output_dir
    mv card*.fasta output_dir

    DB_VERSION=\$(ls $card/card_database_*_all.fasta | sed "s/$card\\/card_database_v\\([0-9].*[0-9]\\).*/\\1/")

    RGI_VERSION=\$(rgi main --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """

stub:
    """
    touch card.fasta
    touch card_all.fasta

    mkdir output_dir
    mv card*.fasta output_dir

    RGI_VERSION=\$(rgi main --version)
    DB_VERSION=stub_version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """
}
