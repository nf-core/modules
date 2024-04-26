#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENRICH } from '../../../../modules/nf-core/genrich/main.nf'
include { GENRICH as GENRICH_CTRL    } from '../../../../modules/nf-core/genrich/main.nf'
include { GENRICH as GENRICH_SE      } from '../../../../modules/nf-core/genrich/main.nf'
include { GENRICH as GENRICH_ALL     } from '../../../../modules/nf-core/genrich/main.nf'
include { GENRICH as GENRICH_ATACSEQ } from '../../../../modules/nf-core/genrich/main.nf'
include { GENRICH as GENRICH_LIST } from '../../../../modules/nf-core/genrich/main.nf'

workflow test_genrich {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  [ ]]
    blacklist = [ ]

    GENRICH ( input, blacklist )
}

workflow test_genrich_ctrl {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    blacklist = [ ]

    GENRICH_CTRL ( input, blacklist )
}

workflow test_genrich_se {
    input     = [ [ id:'test', single_end:true ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    blacklist = [ ]

    GENRICH_SE ( input, blacklist )
}


workflow test_genrich_all_outputs {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    blacklist = [ ]

    GENRICH_ALL ( input, blacklist )
}

workflow test_genrich_blacklist {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  [ ]]
    blacklist = [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]

    GENRICH ( input, blacklist )
}

workflow test_genrich_atacseq {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  [ ]]
    blacklist = [ ]

    GENRICH_ATACSEQ ( input, blacklist )
}

workflow test_genrich_list {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true),
                    file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true)],
                  [ ]]
    blacklist = [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]

    GENRICH_LIST ( input, blacklist )
}

