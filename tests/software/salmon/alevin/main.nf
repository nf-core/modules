#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_ALEVIN } from '../../../../software/salmon/alevin/main.nf' addParams( options: [:] )

workflow test_salmon_alevin {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    SALMON_ALEVIN ( input )
}
