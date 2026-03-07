process BRACKEN_BUILD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f30aa99d8d4f6ff1104f56dbacac95c1dc0905578fb250c80f145b6e80703bd1/data'
        : 'community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01'}"

    input:
    tuple val(meta), path(k2d, stageAs: "bracken-database/"), path(map, stageAs: "bracken-database/"), path(library, stageAs: "bracken-database/library/added/"), path(taxonomy, stageAs: "bracken-database/taxonomy/")

    output:
    tuple val(meta), path("bracken-database/", includeInputs: true), emit: db
    tuple val(meta), path("bracken-database/database*", includeInputs: true), path("bracken-database/*k2d", includeInputs: true), path("bracken-database/*map", includeInputs: true), path("bracken-database/library/added/*", includeInputs: true), path("bracken-database/taxonomy/*", includeInputs: true), emit: db_separated
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
        -d bracken-database/

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
    mkdir bracken-database/
    mkdir -p bracken-database/library/added bracken-database/taxonomy
    touch bracken-database/database100mers.kmer_distrib
    touch bracken-database/database100mers.kraken
    touch bracken-database/database.kraken
    touch bracken-database/{hash,opts,tax}.k2d
    touch bracken-database/seqid2taxid.map
    touch bracken-database/library/added/test.fa
    touch bracken-database/taxonomy/{nodes,names}.dmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """
}
