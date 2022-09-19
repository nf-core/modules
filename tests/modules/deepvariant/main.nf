#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPVARIANT                          } from '../../../modules/deepvariant/main.nf'
include { DEEPVARIANT as DEEPVARIANT_INTERVALS } from '../../../modules/deepvariant/main.nf'

workflow test_deepvariant {

    bam_tuple_ch = Channel.of([ [ id:'test', single_end:false ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                               []
                                ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DEEPVARIANT ( bam_tuple_ch, fasta, fai)
}

workflow test_deepvariant_cram_intervals {

    cram_tuple_ch = Channel.of([[ id:'test', single_end:false ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
                            ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DEEPVARIANT_INTERVALS ( cram_tuple_ch, fasta, fai)
}
