process SUPPA_CLUSTEREVENTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    tuple val(meta), path(dpsi), path(psivec)
    val significance_threshold
    val dpsi_threshold
    val maximum_distance
    val metric
    val separation
    val min_number_events
    val groups
    val clustering_method

    output:
    tuple val(meta), path("*.clustvec"), emit: clustvec
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sig_threshold_arg = significance_threshold ? "--sig-threshold ${significance_threshold}" : ''
    def dpsi_threshold_arg = dpsi_threshold ? "--dpsi-threshold ${dpsi_threshold}" : ''
    def eps_arg = maximum_distance ? "--eps ${maximum_distance}" : ''
    def metric_arg = metric ? ( ["euclidean", "manhattan", "cosine"].contains(metric) ? "--metric ${metric}" : error("Invalid metric: ${metric}. Must be one of: euclidean, manhattan, cosine") ) : ''
    def separation_arg = separation ? "--separation ${separation}" : ''
    def min_pts_arg = min_number_events ? "--min-pts ${min_number_events}" : ''
    def clustering_arg = clustering_method ? ( [ "DBSCAN", "OPTICS" ].contains(clustering_method) ? "--clustering ${clustering_method}" : error("Invalid clustering method: ${clustering_method}. Must be one of: DBSCAN, OPTICS") ) : ''

    clustering_method == 'optics' && !separation ? error("The 'optics' clustering method requires the '--separation' parameter to be set.") : ''

    """
    touch ${prefix}.clustvec

    suppa.py \\
        clusterEvents \\
        --dpsi ${dpsi} \\
        --psivec ${psivec} \\
        ${sig_threshold_arg} \\
        ${dpsi_threshold_arg} \\
        ${eps_arg} \\
        ${metric_arg} \\
        ${separation_arg} \\
        ${min_pts_arg} \\
        --groups ${groups} \\
        ${clustering_arg} \\
        --output ${prefix} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.clustvec
    """
}
