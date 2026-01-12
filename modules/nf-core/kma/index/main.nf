process KMA_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fc6c961562aef21c24b4f2330d9cd7e9bbda162b0d584a5cd5428e0b725e0d6/data':
        'community.wave.seqera.io/library/kma:1.5.0--eb093e0381fb59ea' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("kmaindex"),  emit: index
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${fasta.baseName}"
    def args    = task.ext.args ?: ''
    """
    mkdir kmaindex
    kma \\
        index \\
        -i ${fasta} \\
        -o kmaindex/${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma_index -v 2>&1) | sed 's/^KMA_index-\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir kmaindex

    touch kmaindex/${prefix}.comp.b
    touch kmaindex/${prefix}.length.b
    touch kmaindex/${prefix}.name
    touch kmaindex/${prefix}.seq.b

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma_index -v 2>&1) | sed 's/^KMA_index-\$//')
    END_VERSIONS
    """
}
