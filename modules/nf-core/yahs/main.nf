process YAHS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yahs:1.2.2--h577a1d6_1':
        'biocontainers/yahs:1.2.2--h577a1d6_1' }"

    input:
    tuple val(meta), path(fasta), path(fai), path(hic_map)

    output:
    tuple val(meta), path("*_scaffolds_final.fa") , emit: scaffolds_fasta,  optional: true
    tuple val(meta), path("*_scaffolds_final.agp"), emit: scaffolds_agp  ,  optional: true
    tuple val(meta), path("*.bin")                , emit: binary
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    yahs \\
        -o ${prefix} \\
        ${args} \\
        ${fasta} \\
        ${hic_map}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yahs: \$(yahs --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_scaffolds_final.fa
    touch ${prefix}_scaffolds_final.agp
    touch ${prefix}.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yahs: \$(yahs --version 2>&1)
    END_VERSIONS
    """
}
