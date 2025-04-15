process DSSP_MKDSSP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d47606b41c3a33a667486f98cdbfde7602c8e476bb78eae78d9bf6c86bec5c7b/data':
        'community.wave.seqera.io/library/dssp:4.4.11--d5cd36c6e251360d' }"

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
        --output-format=dssp \\
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
