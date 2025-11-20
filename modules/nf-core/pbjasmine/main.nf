process PBJASMINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbjasmine:2.4.0--h9948957_1':
        'biocontainers/pbjasmine:2.4.0--h9948957_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    jasmine \\
        $args \\
        --num-threads $task.cpus \\
        $bam \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbjasmine: \$(jasmine --version | head -n1 | sed 's/jasmine //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbjasmine: \$(jasmine --version | head -n1 | sed 's/jasmine //')
    END_VERSIONS
    """
}
