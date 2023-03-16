#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_UNIONBEDG                             } from '../../../../../modules/nf-core/bedtools/unionbedg/main.nf'

workflow test_bedtools_unionbedg {
    
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true),
		        file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true) 
		      ]
            ]

    BEDTOOLS_UNIONBEDG ( input )
}
