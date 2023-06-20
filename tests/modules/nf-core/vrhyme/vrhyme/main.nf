#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VRHYME_VRHYME } from '../../../../../modules/nf-core/vrhyme/vrhyme/main.nf'

workflow test_vrhyme_vrhyme {

    bam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    VRHYME_VRHYME ( bam, fasta )
}
