#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADAPTERREMOVAL                            } from '../../../modules/adapterremoval/main.nf'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_COLLAPSE } from '../../../modules/adapterremoval/main.nf'


workflow test_adapterremoval_single_end {
    input = [ [ id:'test', single_end:true, collapse:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]

    ADAPTERREMOVAL ( input, [] )
}

workflow test_adapterremoval_paired_end {
    input = [ [ id:'test', single_end:false, collapse:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    ADAPTERREMOVAL ( input, [] )
}

workflow test_adapterremoval_paired_end_collapse {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    ADAPTERREMOVAL_COLLAPSE ( input, [] )
}

