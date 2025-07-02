process MULTISEQDEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-seurat_r-seuratobject:4c5a804804327d29':
        'community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827' }"

    input:
    tuple val(meta), path(seurat_object), val(assay)

    output:
    tuple val(meta), path("*_params_multiseqdemux.csv"), emit: params
    tuple val(meta), path("*_res_multiseqdemux.csv")   , emit: results
    tuple val(meta), path("*_multiseqdemux.rds")       , emit: rds
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    quantile   = task.ext.assay      ?: 0.7          // The quantile to use for classification
    autoThresh = task.ext.autoThresh ?: true         // Whether to perform automated threshold finding to define the best quantile
    maxiter    = task.ext.maxiter    ?: 5            // Maximum number of iterations if autoThresh is true
    qrangeFrom = task.ext.qrangeFrom ?: 0.1          // A range of possible quantile values to try if autoThresh is true
    qrangeTo   = task.ext.qrangeTo   ?: 0.9          // A range of possible quantile values to try if autoThresh is true
    qrangeBy   = task.ext.qrangeBy   ?: 0.05         // A range of possible quantile values to try if autoThresh is true
    verbose    = task.ext.verbose    ?: true         // Prints the output
    prefix     = task.ext.prefix     ?: "${meta.id}" // Prefix name for the file containing the output of MULTI-Seq Demux object

    template 'MultiSeqDemux.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_params_multiseqdemux.csv
    touch ${prefix}_res_multiseqdemux.csv
    touch ${prefix}_multiseqdemux.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
    END_VERSIONS
    """
}
