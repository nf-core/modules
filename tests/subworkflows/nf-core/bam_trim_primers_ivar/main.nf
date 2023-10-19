#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_TRIM_PRIMERS_IVAR } from '../../../../subworkflows/nf-core/bam_trim_primers_ivar/main.nf'

workflow test_bam_trim_primers_ivar {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    bed_file   = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BAM_TRIM_PRIMERS_IVAR ( input, bed_file, fasta )
}
