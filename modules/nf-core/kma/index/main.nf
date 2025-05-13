process KMA_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.15--h577a1d6_1' :
        'biocontainers/kma:1.4.15--h577a1d6_1' }"

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
