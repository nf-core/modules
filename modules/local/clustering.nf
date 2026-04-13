process CLUSTERING {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.python_container}"
    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    tuple val(meta), path(pca_scores), path(pca_info)
    path clustering_script

    output:
    tuple val(meta), path('clusters.csv'),         emit: clusters
    tuple val(meta), path('clustering_info.json'), emit: clustering_info

    script:
    def algorithm = params.algorithm ?: 'kmeans'
    def n_init    = params.n_init ?: 100
    def init_meth = params.init_method ?: 'random'

    def extra_args = ''
    if( algorithm == 'kmeans' ) {
        def k = params.n_clusters ?: 3
        extra_args = "--k ${k} --n_init ${n_init} --init-method ${init_meth}"
    }
    else if( algorithm == 'dbscan' ) {
        def eps  = params.dbscan_eps
        def mins = params.dbscan_min_samples
        extra_args = "--dbscan-eps ${eps} --dbscan-min-samples ${mins}"
    }

    """
    python3 ${clustering_script} \\
        --features ${pca_scores} \\
        --algorithm ${algorithm} \\
        ${extra_args} \\
        --out-clusters clusters.csv \\
        --out-info clustering_info.json
    """
}