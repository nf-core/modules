process RGI_CARDANNOTATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f452c8e124ee58ab6b26442d15401c57d471cb753f53921570dc484df4e7620/data'
        : 'community.wave.seqera.io/library/rgi_kma:e905ecb8305e2609' }"

    input:
    path card

    output:
    path ("card_database_processed"), emit: db
    env 'RGI_VERSION', emit: tool_version
    env 'DB_VERSION', emit: db_version
    tuple val("${task.process}"), val('rgi'), eval("rgi main --version"),  emit: versions_rgi, topic: versions
    tuple val("${task.process}"), val('rgi-database'), eval("echo \$DB_VERSION"),  emit: versions_db , topic: versions
    tuple val("${task.process}"), val('kma'), eval("kma -v"),  emit: versions_kma, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    rgi card_annotation \\
        -i ${card}/card.json \\
        ${args}

    DB_VERSION=\$(ls card_database_*_all.fasta | sed "s/card_database_v\\([0-9].*[0-9]\\).*/\\1/")

    mkdir card_database_processed
    mv card*.fasta card_database_processed
    cp ${card}/* card_database_processed

    RGI_VERSION=\$(rgi main --version)
    """

    stub:
    """
    touch card.fasta
    touch card_all.fasta

    mkdir card_database_processed
    mv card*.fasta card_database_processed

    RGI_VERSION=\$(rgi main --version)
    DB_VERSION=stub_version
    """
}
