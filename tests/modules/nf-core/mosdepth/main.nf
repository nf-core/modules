#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOSDEPTH                       } from '../../../../modules/nf-core/mosdepth/main.nf'
include { MOSDEPTH as MOSDEPTH_FAIL      } from '../../../../modules/nf-core/mosdepth/main.nf'
include { MOSDEPTH as MOSDEPTH_WINDOW    } from '../../../../modules/nf-core/mosdepth/main.nf'
include { MOSDEPTH as MOSDEPTH_THRESHOLD } from '../../../../modules/nf-core/mosdepth/main.nf'
include { MOSDEPTH as MOSDEPTH_QUANTIZED } from '../../../../modules/nf-core/mosdepth/main.nf'

workflow test_mosdepth {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    MOSDEPTH ( input, [[:],[]], [[:],[]] )
}

workflow test_mosdepth_bed {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    bed = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    MOSDEPTH ( input, bed, [[:],[]] )
}

workflow test_mosdepth_cram {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]
    fasta = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MOSDEPTH ( input, [[:],[]], fasta )
}

workflow test_mosdepth_cram_bed {
    input  = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]
    bed = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]
    fasta = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    MOSDEPTH ( input, bed, fasta )
}

workflow test_mosdepth_window {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    bed = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    MOSDEPTH_WINDOW ( input, [[:],[]], [[:],[]] )
}

workflow test_mosdepth_quantized {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    MOSDEPTH_QUANTIZED ( input, [[:],[]], [[:],[]] )
}

workflow test_mosdepth_thresholds {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    bed = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]
    MOSDEPTH_THRESHOLD ( input, bed, [[:],[]] )
}

workflow test_mosdepth_fail {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    bed = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]
    MOSDEPTH_FAIL ( input, bed, [[:],[]] )
}
