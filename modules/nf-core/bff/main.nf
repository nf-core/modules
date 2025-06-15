process BFF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-seurat:3.0.2--r36h0357c0b_0':
        'biocontainers/r-seurat:3.0.2--r36h0357c0b_0' }"

    input:
    tuple val(meta), path(hto_matrix), val(methodsForConsensus), val(preprocess_bff)

    output:
    tuple val(meta), path("*_assignment_bff.csv")       , emit: results
    tuple val(meta), path("*_params_bff.csv")           , emit: params
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"


    bff = "True"
    rna_matrix_bff = "raw"
    hto_matrix_bff = "raw"
    assignmentOutBff = "bff"
    bff_preprocess = "True"
    methods = "combined_bff"
    methodsForConsensus = 'NULL'
    cellbarcodeWhitelist = 'NULL'
    metricsFile = 'metrics_bff.csv'
    doTSNE = "True"
    doHeatmap = "True"
    perCellSaturation = 'NULL'
    majorityConsensusThreshold = 'NULL'
    chemistry = "10xV3"
    callerDisagreementThreshold = 'NULL'
    preprocess_bff = "FALSE"
    barcodeWhitelist = 'NULL'


            def run_preprocess = preprocess_bff != 'False' ? " --preprocess_bff" : ''
        """
        mkdir bff_${sampleId}
        bff.R --fileHto hto_data --methods $methods --methodsForConsensus $methodsForConsensus \
        --cellbarcodeWhitelist $cellbarcodeWhitelist --metricsFile bff_${sampleId}_$metricsFile \
        --doTSNE $doTSNE --doHeatmap $doHeatmap --perCellSaturation $perCellSaturation --majorityConsensusThreshold $majorityConsensusThreshold \
        --chemistry $chemistry --callerDisagreementThreshold $callerDisagreementThreshold --outputdir bff_${sampleId} \
        --assignmentOutBff $assignmentOutBff ${run_preprocess} --barcodeWhitelist $barcodeWhitelist
        """




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
