process AUTHENTICT_DEAM2CONT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/authentict:1.0.1--py311h9f5acd7_0':
        'biocontainers/authentict:1.0.1--py311h9f5acd7_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(config)
    tuple val(meta3), path(positions)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config_file = config ? "-c ${config}" : ""
    def positions_file = positions ? "-p ${positions}" : ""

    """
    samtools view $args $bam | AuthentiCT \\
        deam2cont \\
        $args2 \\
        $config_file \\
        $positions_file \\
        - \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        authentict: \$(echo \$(AuthentiCT --version 2>&1) )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub :
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        authentict: \$(echo \$(AuthentiCT --version 2>&1) )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
