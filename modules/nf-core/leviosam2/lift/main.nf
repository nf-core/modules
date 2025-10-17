process LEVIOSAM2_LIFT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leviosam2:0.4.2--h4ac6f70_0':
        'biocontainers/leviosam2:0.4.2--h4ac6f70_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta_ref), path(clft)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    leviosam2 \\
        lift \\
        -t ${task.cpus} \\
        -C ${clft} \\
        -a ${input} \\
        -p ${prefix} \\
        -O bam \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leviosam2: \$(leviosam2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leviosam2: \$(leviosam2 --version)
    END_VERSIONS
    """
}
