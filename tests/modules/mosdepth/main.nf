#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOSDEPTH } from '../../../modules/mosdepth/main.nf'

workflow test_mosdepth {
    input  = [ [ id:'test', single_end:true ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
             ]

    MOSDEPTH ( input, [], [] )
}


workflow test_mosdepth_window {
    input  = [ [ id:'test', single_end:true ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
             ]
    window = 100

    MOSDEPTH ( input, [], window )
}


workflow test_mosdepth_bed {
    input  = [ [ id:'test', single_end:true ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
             ]
    bed  = [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]

    MOSDEPTH ( input, bed, [] )
}
