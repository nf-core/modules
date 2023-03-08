#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GOLEFT_INDEXSPLIT } from '../../../../../modules/nf-core/goleft/indexsplit/main.nf'

workflow test_goleft_indexsplit {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fai = [
        [ id:'sarscov2'], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    ]

    GOLEFT_INDEXSPLIT ( input, fai, 10 )
}
