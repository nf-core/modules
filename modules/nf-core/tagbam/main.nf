process TAGBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c4d2c045368ba77115020250c39205495f7e42f9d49e4f98add77280757cc2c/data':
        'community.wave.seqera.io/library/tagbam:0.1.0--a58e5b41b01b9027' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    tagbam \\
        $args \\
        --threads $task.cpus \\
        --input $bam \\
        --output-file ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tagbam: \$(tagbam --version | sed 's/tagbam //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tagbam: \$(tagbam --version | sed 's/tagbam //')
    END_VERSIONS
    """
}
