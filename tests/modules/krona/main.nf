#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA } from '../../../modules/krona/main.nf' addParams( options: [:] )

workflow test_krona {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    KRONA ( input )
}
