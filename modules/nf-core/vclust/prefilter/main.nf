process VCLUST_PREFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vclust \\
        prefilter \\
        $args \\
        -t $task.cpus \\
        -i ${fasta} \\
        -o ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version)
    END_VERSIONS
    """
}
