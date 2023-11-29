#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_ZIPPERBAMS } from '../../../../../modules/nf-core/fgbio/zipperbams/main.nf'

workflow test_fgbio_zipperbams {

    mbam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_mapped_bam'], checkIfExists: true)
    ]

    ubam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_unmapped_bam'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'],                          checkIfExists: true)
    dict  = file(params.test_data['homo_sapiens']['genome']['genome_dict'],                           checkIfExists: true)


    FGBIO_ZIPPERBAMS ( ubam, mbam, fasta, dict )
}
