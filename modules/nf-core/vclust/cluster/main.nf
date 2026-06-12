process VCLUST_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'quay.io/biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    input:
    tuple val(meta), path(tsv)
    tuple val(meta2), path(ids)
    val metric
    val tani
    val gani
    val ani

    output:
    tuple val(meta), path("*.tsv"), emit: clusters
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('vclust'), eval("vclust --version | sed 's/v//'"), topic: versions, emit: versions_vclust


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def tani_command = tani ? "--tani ${tani} " : ''
    def gani_command = gani ? "--gani ${gani} " : ''
    def ani_command = ani ? "--ani ${ani} " : ''
    def metric_command = metric ? "--metric ${metric} " : ''
    def tgani_command = tani_command + gani_command + ani_command ?: ''

    if ( !metric_command ) error "ERROR: The metric must be specified."
    if ( !tgani_command ) error "ERROR: At least one of the ani thresholds must be specified: --tani, --gani, --ani."
    if ( !tgani_command.contains("--${metric}")) error "ERROR: The metric '${metric}' must have an associated threshold."

    """
    vclust \\
        cluster \\
        ${args} \\
        ${metric_command} \\
        ${tgani_command} \\
        -i ${tsv} \\
        --ids ${ids} \\
        -o ${prefix}.cluster.tsv 2>&1 | tee ${prefix}.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tani_command = tani ? "--tani ${tani}" : ''
    def gani_command = gani ? "--gani ${gani}" : ''
    def ani_command = ani ? "--ani ${ani}" : ''
    def metric_command = metric ? "--metric ${metric}" : ''
    def tgani_command = tani_command + gani_command + ani_command ?: ''

    if ( !metric_command ) error "ERROR: The metric must be specified."
    if ( !tgani_command ) error "ERROR: At least one of the ani thresholds must be specified: --tani, --gani, --ani."
    if ( !tgani_command.contains("${metric}")) error "ERROR: The metric '${metric}' must have an associated threshold."

    """
    touch ${prefix}.clusters.tsv
    echo "${metric_command} ${tgani_command}" > ${prefix}.log
    """
}
