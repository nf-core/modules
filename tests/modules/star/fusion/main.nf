#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_FUSION } from '../../../../modules/star/fusion/main.nf' addParams( options: [:] )

workflow test_star_fusion {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    STAR_FUSION ( input )
}
