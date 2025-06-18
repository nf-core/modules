process HASHEDDROPS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-dropletutils_r-seurat:3cbdf18d48cd0cfa':
        'community.wave.seqera.io/library/bioconductor-dropletutils_r-seurat:e1dff3a0fb7c5920' }"

    input:
    tuple val(meta), path(hto_matrix), val(runEmptyDrops), path(rna_matrix)

    output:
    tuple val(meta), path("*_emptyDrops.png")          , emit: empty_drops_plot
    tuple val(meta), path("*_emptyDrops.csv")          , emit: empty_drops_csv
    tuple val(meta), path("*_emptyDrops.rds")          , emit: empty_drops_rds
    tuple val(meta), path("*_results_hasheddrops.csv") , emit: results
    tuple val(meta), path("*_hasheddrops.rds")         , emit: rds
    tuple val(meta), path("*_plot_hasheddrops.png")    , emit: plot
    tuple val(meta), path("*_params_hasheddrops.csv")  , emit: params
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // emptyDrops Parameters
    lower           = task.ext.lower           ?: "100"        // A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
    niters          = task.ext.niters          ?: "10000"      // An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
    testAmbient     = task.ext.testAmbient     ?: "TRUE"       // A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower.
    round           = task.ext.round           ?: "TRUE"       // Logical scalar indicating whether to check for non-integer values in m and, if present, round them for ambient profile estimation.
    byRank          = task.ext.byRank          ?: "NULL"       // An integer scalar parametrizing an alternative method for identifying assumed empty droplets. If set, this is used to redefine lower and any specified value for lower is ignored.
    isCellFDR       = task.ext.isCellFDR       ?: "0.01"       // Threshold to filter the cells.
    gene_col        = task.ext.gene_col        ?: "2"          // Specify which column of genes.tsv or features.tsv to use for gene names; default is 2.

    // hashedDrops Parameters
    ignore          = task.ext.ignore          ?: "NULL"       // A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored.
    alpha           = task.ext.alpha           ?: "NULL"       // A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme.
    ambient         = task.ext.ambient         ?: "TRUE"       // Whether to use the relative abundance of each HTO in the ambient solution from emptyDrops, set TRUE only when test_ambient is TRUE.
    minProp         = task.ext.minProp         ?: "0.05"       // Numeric scalar to be used to infer the ambient profile when ambient=NULL.
    pseudoCount     = task.ext.pseudoCount     ?: "5"          // A numeric scalar specifying the minimum pseudo-count when computing logfold changes.
    constantAmbient = task.ext.constantAmbient ?: "FALSE"      // Logical scalar indicating whether a constant level of ambient contamination should be used to estimate LogFC2 for all cells.
    doubletNmads    = task.ext.doubletNmads    ?: "3"          // A numeric scalar specifying the number of median absolute deviations (MADs) to use to identify doublets.
    doubletMin      = task.ext.doubletMin      ?: "2"          // A numeric scalar specifying the minimum threshold on the log-fold change to use to identify doublets.
    doubletMixture  = task.ext.doubletMixture  ?: "FALSE"      // Logical scalar indicating whether to use a 2-component mixture model to identify doublets.
    confidentNmads  = task.ext.confidentNmads  ?: "3"          // A numeric scalar specifying the number of MADs to use to identify confidently assigned singlets.
    confidentMin    = task.ext.confidentMin    ?: "2"          // A numeric scalar specifying the minimum threshold on the log-fold change to use to identify singlets.
    combinations    = task.ext.combinations    ?: "NULL"       // An integer matrix specifying valid combinations of HTOs. Each row corresponds to a single sample and specifies the indices of rows in x corresponding to the HTOs used to label that sample.

    // others
    prefix          = task.ext.prefix          ?: "${meta.id}" // Prefix name for output files.

    template 'HashedDrops.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_emptyDrops.png
    touch ${prefix}_emptyDrops.csv
    touch ${prefix}_emptyDrops.rds
    touch ${prefix}_results_hasheddrops.csv
    touch ${prefix}_hasheddrops.rds
    touch ${prefix}_plot_hasheddrops.png
    touch ${prefix}_params_hasheddrops.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        cdropletutils: \$(Rscript -e "library(DropletUtils); cat(as.character(packageVersion('DropletUtils')))")
    END_VERSIONS
    """
}

