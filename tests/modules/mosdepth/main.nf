#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOSDEPTH                  } from '../../../modules/mosdepth/main.nf'
include { MOSDEPTH as MOSDEPTH_FAIL } from '../../../modules/mosdepth/main.nf'

workflow test_mosdepth {
    input  = [
                [ id:'test', single_end:true ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
            ]

    MOSDEPTH ( input, [], [] )
}

workflow test_mosdepth_bed {
    input  = [
                [ id:'test', single_end:true ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
            ]
    bed  = [ file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true) ]

    MOSDEPTH ( input, bed, [] )
}

workflow test_mosdepth_cram {
    input  = [
                [ id:'test', single_end:true ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true) ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true) ]
            ]
    fasta = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]

    MOSDEPTH ( input, [], fasta )
}

workflow test_mosdepth_cram_bed {
    input  = [
                [ id:'test', single_end:true ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true) ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true) ]
            ]
    bed  = [ file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true) ]
    fasta = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]

    MOSDEPTH ( input, bed, fasta )
}

workflow test_mosdepth_fail {
    input  = [
                [ id:'test', single_end:true ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
                [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
            ]
    bed  = [ file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true) ]

    MOSDEPTH_FAIL ( input, bed, [] )
}
