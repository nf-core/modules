process BFF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-seurat:3.0.2--r36h0357c0b_0':
        'biocontainers/r-seurat:3.0.2--r36h0357c0b_0' }"

    input:
    tuple val(meta), path(hto_matrix), val(methods), val(preprocessing)

    output:
    tuple val(meta), path("*_assignment_bff.csv")       , emit: results
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

    template bff.R

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-propr: \$(Rscript -e "library('propr'); cat(as.character(packageVersion('propr')))")
    END_VERSIONS
    """
}
