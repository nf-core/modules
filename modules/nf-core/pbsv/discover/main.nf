process PBSV_DISCOVER {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0a49662ef2db75cd0afdf7f90840bc82f8f96315af1ffb2e2391e55a0ff1c861/data':
        'community.wave.seqera.io/library/pbsv:2.11.0--c85e7f17330a07c9' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.svsig.gz"), emit: svsig
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbsv \\
        discover \\
        $args \\
        ${bam} \\
        ${prefix}.svsig.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsv: \$(pbsv --version |& sed '1!d ; s/pbsv //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.svsig.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsv: \$(pbsv --version |& sed '1!d ; s/pbsv //')
    END_VERSIONS
    """
}
