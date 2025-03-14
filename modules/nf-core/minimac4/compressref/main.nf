process MINIMAC4_COMPRESSREF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimac4:4.1.6--hcb620b3_1':
        'biocontainers/minimac4:4.1.6--hcb620b3_1' }"

    input:
    tuple val(meta), path(ref), path(ref_index) // Reference index is autodetected from reference file name

    output:
    tuple val(meta), path("*.msav"), emit: msav
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minimac4 \\
        --compress-reference $ref\\
        $args \\
        --threads $task.cpus \\
        -o ${prefix}.msav \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: \$(minimac4 --version |& sed '1!d ; s/minimac v//')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.msav

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: \$(minimac4 --version |& sed '1!d ; s/minimac v//')
    END_VERSIONS
    """
}
