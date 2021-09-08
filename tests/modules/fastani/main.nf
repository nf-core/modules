#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTANI } from '../../../modules/fastani/main.nf' addParams( options: [:] )

workflow test_fastani {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    FASTANI ( input )
}
