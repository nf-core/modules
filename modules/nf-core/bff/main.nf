process BFF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-cellhashr_r-seurat:1c94360d8ed188c4':
        'community.wave.seqera.io/library/bioconductor-cellhashr_r-seurat:25c4bc76749af5ac' }"

    input:
    tuple val(meta), path(hto_matrix), val(methods), val(preprocessing)

    output:
    tuple val(meta), path("*_assignment_bff.csv"), emit: assignment
    tuple val(meta), path("*_metrics_bff.csv")   , emit: metrics
    tuple val(meta), path("*_params_bff.csv")    , emit: params
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    methodsForConsensus         = task.ext.methodsForConsensus         ?: null       // By default, a consensus call will be generated using all methods, null, RAW OR CLUSTER
    cellbarcodeWhitelist        = task.ext.cellbarcodeWhitelist        ?: null       // A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix
    prefix                      = task.ext.prefix                      ?: "${meta.id}" // Prefix name for output files
    doTSNE                      = task.ext.doTSNE                      ?: false       // If true, tSNE will be run on the resulting hashing calls after each caller
    barcodeWhitelist            = task.ext.barcodeWhitelist            ?: null       // A vector of barcode names to retain, used for preprocessing step
    doHeatmap                   = task.ext.doHeatmap                   ?: true       // If true, Seurat::HTOHeatmap will be run on the results of each caller
    perCellSaturation           = task.ext.perCellSaturation           ?: null       // An optional dataframe with the columns cellbarcode and saturation
    majorityConsensusThreshold  = task.ext.majorityConsensusThreshold  ?: null       // This applies to calculating a consensus call when multiple algorithms are used
    chemistry                   = task.ext.chemistry                   ?: "10xV3"      // This string is passed to EstimateMultipletRate. Should be either 10xV2 or 10xV3
    callerDisagreementThreshold = task.ext.callerDisagreementThreshold ?: null       // If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells

    template 'bff.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_assignment_bff.csv
    touch ${prefix}_metrics_bff.csv
    touch ${prefix}_params_bff.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        cellhashR: \$(Rscript -e "library(cellhashR); cat(as.character(packageVersion('cellhashR')))")
    END_VERSIONS
    """
}
