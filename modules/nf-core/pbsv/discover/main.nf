process PBSV_DISCOVER {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pbsv:2.11.0--64a2437cb4119823':
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
