#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHROMAP } from '../../../software/chromap/main.nf' addParams( options: [:] )

workflow test_chromap {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    CHROMAP ( input )
}
