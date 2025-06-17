process BFF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-cellhashr_r-seurat:1c94360d8ed188c4':
        'community.wave.seqera.io/library/bioconductor-cellhashr_r-seurat:25c4bc76749af5ac' }"

    input:
    tuple val(meta), path(hto_matrix), val(methods), val(preprocessing)

    output:
    tuple val(meta), path("*_assignment_bff.csv")       , emit: assignment
    tuple val(meta), path("*_metrics_bff.csv")          , emit: metrics
    tuple val(meta), path("*_params_bff.csv")           , emit: params
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    methodsForConsensus = task.ext.methodsForConsensus ?: "NULL" // NULL, RAW OR CLUSTER
    cellbarcodeWhitelist = task.ext.cellbarcodeWhitelist ?: "NULL"
    prefix = task.ext.prefix ?: "${meta.id}"

    doTSNE = task.ext.doTSNE ?: "TRUE"
    barcodeWhitelist = task.ext.barcodeWhitelist ?: "NULL"
    doHeatmap = task.ext.doHeatmap ?: "TRUE"
    perCellSaturation = task.ext.perCellSaturation ?: "NULL"
    majorityConsensusThreshold = task.ext.majorityConsensusThreshold ?: "NULL"
    chemistry = task.ext.chemistry ?: "10xV3"
    callerDisagreementThreshold = task.ext.callerDisagreementThreshold ?: "NULL"

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
