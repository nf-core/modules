#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATAQV_ATAQV } from '../../../../modules/ataqv/ataqv/main.nf'
include { ATAQV_ATAQV as ATAQV_ATAQV_PROBLEM_READS} from '../../../../modules/ataqv/ataqv/main.nf'

workflow test_ataqv_ataqv {
    
    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        [],
        []
    ]
    
    ATAQV_ATAQV ( input, 'human', [], [], [] )
}

workflow test_ataqv_ataqv_problem_reads {
    
    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        [],
        []
    ]
    
    ATAQV_ATAQV_PROBLEM_READS ( input, 'human', [], [], [] )
}

workflow test_ataqv_ataqv_peak {
    
    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        [],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]
    
    ATAQV_ATAQV ( input, 'human', [], [], [] )
}

workflow test_ataqv_ataqv_tss {
    
    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        []
    ] 
    tss_file = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
   
    ATAQV_ATAQV ( input, 'human', tss_file, [], [] )
}

workflow test_ataqv_ataqv_excluded_regs {
    
    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        []
    ] 
    tss_file = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    excl_regs_file = file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
   
    ATAQV_ATAQV ( input, 'human', tss_file, excl_regs_file, [] )
}
