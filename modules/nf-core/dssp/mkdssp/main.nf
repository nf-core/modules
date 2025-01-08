process DSSP_MKDSSP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/62e00156fdddc2b7933c1221e69715a2393813fb24f68a0cdc0a4fc949daad2a/data':
        'community.wave.seqera.io/library/dssp:b16ae164ced98723' }"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("*.{dssp,mmcif}"), emit: dssp
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdssp \\
        $args \\
        ${pdb} \\
        ${prefix}.dssp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dssp: \$(dssp --version |& sed '1!d ; s/dssp //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dssp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dssp: \$(dssp --version |& sed '1!d ; s/dssp //')
    END_VERSIONS
    """
}
