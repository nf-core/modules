process YAK_COUNT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yak:0.1--he4a0461_4':
        'biocontainers/yak:0.1--he4a0461_4' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.yak"), emit: yak
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input_command = meta.single_end ? "${fastq}" : "<(zcat ${fastq}) <(zcat ${fastq})"
    """
    yak \\
        count \\
        $args \\
        -t${task.cpus} \\
        -o ${prefix}.yak \\
        $input_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yak: \$(yak version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yak

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yak: \$(yak version)
    END_VERSIONS
    """
}
