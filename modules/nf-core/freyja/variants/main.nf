process FREYJA_VARIANTS {
    tag "$meta.id"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e02db57d7bca8a8240a52b1c5e1d8dc4e1ab13e8106ff31d812be86db2be1850/data':
        'community.wave.seqera.io/library/freyja:1.5.3--3b69f13ac540c40e' }"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.variants.tsv"), path("*.depth.tsv"), emit: variants
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    freyja \\
        variants \\
        $args \\
        --ref $fasta \\
        --variants ${prefix}.variants.tsv \\
        --depths ${prefix}.depth.tsv \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variants.tsv
    touch ${prefix}.depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}
