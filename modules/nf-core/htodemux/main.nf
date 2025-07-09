process HTODEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-seurat_r-seuratobject:4c5a804804327d29':
        'community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827' }"

    input:
    tuple val(meta), path(seurat_object), val(assay)

    output:
    tuple val(meta), path("*_params_htodemux.csv")        , emit: params
    tuple val(meta), path("*_assignment_htodemux.csv")    , emit: assignment
    tuple val(meta), path("*_classification_htodemux.csv"), emit: classification
    tuple val(meta), path("*_htodemux.rds")               , emit: rds
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    quantile = task.ext.quantile ?: 0.99         // The quantile of inferred 'negative' distribution for each hashtag, over which the cell is considered positive
    init     = task.ext.init     ?: null         // Initial number of clusters for hashtags
    nstarts  = task.ext.nstarts  ?: 100          // nstarts value for k-means clustering
    kfunc    = task.ext.kfunc    ?: "clara"      // Clustering function for initial hashtag grouping. Default is clara for fast k-medoids clustering on large applications, also support kmeans for kmeans clustering
    nsamples = task.ext.nsamples ?: 100          // Number of samples to be drawn from the dataset used for clustering, for kfunc = clara
    seed     = task.ext.seed     ?: 42           // Sets the random seed
    verbose  = task.ext.verbose  ?: true         // Verbose output
    prefix   = task.ext.prefix   ?: "${meta.id}" // Prefix name for the output files

    template 'HTODemux.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_params_htodemux.csv
    touch ${prefix}_assignment_htodemux.csv
    touch ${prefix}_classification_htodemux.csv
    touch ${prefix}_htodemux.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
