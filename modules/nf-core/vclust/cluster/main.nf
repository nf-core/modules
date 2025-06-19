process VCLUST_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    input:
    tuple val(meta), path(tsv)
    tuple val(meta2), path(ids)
    val tani
    val gani
    val ani

    output:
    tuple val(meta), path("*.tsv"), emit: clusters
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tani_command = tani ? "--tani ${tani}" : ''
    def gani_command = gani ? "--gani ${gani}" : ''
    def ani_command = ani ? "--ani ${ani}" : ''

    def metric_command = tani_command + gani_command + ani_command ?: '--ani 0.95'

    """
    vclust \\
        cluster \\
        $args \\
        ${metric_command} \\
        -i ${tsv} \\
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
