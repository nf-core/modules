#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOSDEPTH } from '../../../software/mosdepth/main.nf' addParams( options: [:] )

workflow test_mosdepth {
    input  = [ [ id:'test', single_end:true ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ] 
             ]
    dummy  = file("dummy_file.txt")
    window = 100

    MOSDEPTH ( input, dummy, window )
}
