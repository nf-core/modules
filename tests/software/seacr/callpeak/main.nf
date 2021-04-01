#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEACR_CALLPEAK } from '../../../../software/seacr/callpeak/main.nf' addParams( options: [ args:'norm stringent' ] )

workflow test_seacr_callpeak {
    input = [ [ id:'test_1'],
              file("${launchDir}/tests/data/generic/bedgraph/K27me3_1_to_chr20.bedgraph", checkIfExists: true),
              file("${launchDir}/tests/data/generic/bedgraph/IgG_1_to_chr20.bedgraph", checkIfExists: true) 
            ]

    SEACR_CALLPEAK ( input )
}
