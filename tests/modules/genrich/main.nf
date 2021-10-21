#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENRICH } from '../../../modules/genrich/main.nf' addParams( control_bam: false, pvalues: false, pileup:false, bed:false, blacklist_bed:false, save_duplicates:false, options: ["args": "-p 0.1"] )
include { GENRICH as GENRICH_BLACKLIST   } from '../../../modules/genrich/main.nf' addParams( control_bam: false, pvalues: false, pileup:false, bed:false, blacklist_bed:true, save_duplicates:false, options: ["args": "-p 0.1"] )
include { GENRICH as GENRICH_ALL_OUTPUTS } from '../../../modules/genrich/main.nf' addParams( control_bam: false, pvalues: true, pileup:true, bed:true, blacklist_bed:false, save_duplicates:true, options: ["args": "-r -p 0.1"] )
include { GENRICH as GENRICH_ATACSEQ     } from '../../../modules/genrich/main.nf' addParams( control_bam: false, pvalues: false, pileup:false, bed:false, blacklist_bed:false, save_duplicates:false, options: ["args": "-j -p 0.1"] )

workflow test_genrich {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ ]
    blacklist = [ ]

    GENRICH ( input, control, blacklist )
}

workflow test_genrich_ctrl {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]
    blacklist = [ ]

    GENRICH ( input, control, blacklist )
}

workflow test_genrich_all_outputs {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]
    blacklist = [ ]

    GENRICH_ALL_OUTPUTS ( input, control, blacklist )
}

workflow test_genrich_atacseq {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]
    blacklist = [ ]

    GENRICH_ATACSEQ ( input, control, blacklist )
}
