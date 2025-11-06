process BRACKEN_BUILD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f30aa99d8d4f6ff1104f56dbacac95c1dc0905578fb250c80f145b6e80703bd1/data'
        : 'community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01'}"

    input:
    tuple val(meta), path(k2d, stageAs: "kraken2db_forbuilding/"), path(map, stageAs: "kraken2db_forbuilding/"), path(library, stageAs: "kraken2db_forbuilding/library/added/"), path(taxonomy, stageAs: "kraken2db_forbuilding/taxonomy/")

    output:
    tuple val(meta), path("${prefix}/", includeInputs: true), emit: db
    tuple val(meta), path("${prefix}/database*", includeInputs: true), emit: bracken_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bracken-build \\
        ${args} \\
        -t ${task.cpus} \\
        -d kraken2db_forbuilding/

    mv kraken2db_forbuilding/ ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    mkdir kraken2db_forbuilding/
    touch kraken2db_forbuilding/database100mers.kmer_distrib
    touch kraken2db_forbuilding/database100mers.kraken
    touch kraken2db_forbuilding/database.kraken

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """
}
