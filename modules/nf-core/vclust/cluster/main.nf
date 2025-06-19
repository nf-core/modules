process VCLUST_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    input:
    tuple val(meta), path(ani)
    tuple val(meta2), path(ids)

    output:
    tuple val(meta), path("*.tsv"), emit: clusters
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vclust \\
        cluster \\
        $args \\
        -i ${ani} \\
        --ids ${ids} \\
        -o ${prefix}.cluster.tsv 2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.clusters.tsv
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version)
    END_VERSIONS
    """
}
