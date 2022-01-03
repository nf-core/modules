#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ALLELECOUNTER } from '../../../modules/allelecounter/main.nf'

workflow test_allelecounter_bam {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    positions = [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]

    ALLELECOUNTER ( input, positions, [] )
}


workflow test_allelecounter_cram {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
            ]
    positions = [ file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true) ]
    fasta = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]

    ALLELECOUNTER ( input, positions, fasta )
}
