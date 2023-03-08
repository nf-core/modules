#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREYJA_VARIANTS } from '../../../../../modules/nf-core/freyja/variants/main.nf'

workflow test_freyja_variants {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    FREYJA_VARIANTS ( input, fasta)
}
