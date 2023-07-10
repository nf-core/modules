#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_UNIONBEDG                             } from '../../../../../modules/nf-core/bedtools/unionbedg/main.nf'
include { BEDTOOLS_UNIONBEDG as BEDTOOLS_UNIONBEDG_EMPTY } from '../../../../../modules/nf-core/bedtools/unionbedg/main.nf'

workflow test_bedtools_unionbedg {
    
    input = [ [ id:'test_output' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true),
		        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) 
		      ]
            ]

    BEDTOOLS_UNIONBEDG ( input, [[:],[]] )
}

workflow test_bedtools_unionbedg_empty {
    
    input = [ [ id:'test_output' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true),
		        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) 
		      ]
            ]
	sizes = [ [ id: 'genome' ],
              file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true) ]

BEDTOOLS_UNIONBEDG_EMPTY ( input, sizes )
}
