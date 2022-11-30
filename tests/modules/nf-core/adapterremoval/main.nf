#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADAPTERREMOVAL                            } from '../../../../modules/nf-core/adapterremoval/main.nf'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_COLLAPSE } from '../../../../modules/nf-core/adapterremoval/main.nf'


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

workflow test_adapterremoval_paired_end_adapterlist {
    input = [ [ id:'test', single_end:false, collapse:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    adapterlist = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/adapterremoval/adapterremoval_adapterlist.txt", checkIfExists: true)

    ADAPTERREMOVAL ( input, adapterlist )
}


