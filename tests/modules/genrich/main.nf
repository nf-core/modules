#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENRICH } from '../../../modules/genrich/main.nf'
include { GENRICH as GENRICH_CTRL    } from '../../../modules/genrich/main.nf'
include { GENRICH as GENRICH_ALL     } from '../../../modules/genrich/main.nf'
include { GENRICH as GENRICH_ATACSEQ } from '../../../modules/genrich/main.nf'

workflow test_genrich {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ ]
    blacklist = [ ]

    save_pvalues    = false
    save_pileup     = false
    save_bed        = false
    save_duplicates = false

    GENRICH ( input, control, blacklist, save_pvalues, save_pileup, save_bed, save_duplicates )
}

workflow test_genrich_ctrl {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]
    blacklist = [ ]

    save_pvalues    = false
    save_pileup     = false
    save_bed        = false
    save_duplicates = false

    GENRICH_CTRL ( input, control, blacklist, save_pvalues, save_pileup, save_bed, save_duplicates )
}

workflow test_genrich_all_outputs {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]
    blacklist = [ ]

    save_pvalues    = true
    save_pileup     = true
    save_bed        = true
    save_duplicates = true

    GENRICH_ALL ( input, control, blacklist, save_pvalues, save_pileup, save_bed, save_duplicates )
}

workflow test_genrich_blacklist {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ ]
    blacklist = [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]

    save_pvalues    = false
    save_pileup     = false
    save_bed        = false
    save_duplicates = false

    GENRICH ( input, control, blacklist, save_pvalues, save_pileup, save_bed, save_duplicates )
}

workflow test_genrich_atacseq {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ]]
    control   = [ ]
    blacklist = [ ]

    save_pvalues    = false
    save_pileup     = false
    save_bed        = false
    save_duplicates = false

    GENRICH_ATACSEQ ( input, control, blacklist, save_pvalues, save_pileup, save_bed, save_duplicates )
}

