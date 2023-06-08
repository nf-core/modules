#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_GREP as SEQKIT_GREP_FILE    } from '../../../../../modules/nf-core/seqkit/grep/main.nf'
include { SEQKIT_GREP as SEQKIT_GREP_STRING  } from '../../../../../modules/nf-core/seqkit/grep/main.nf'

workflow test_seqkit_grep_file {

    sequence = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    pattern = [
        file(params.test_data['homo_sapiens']['genome']['genome_header'], checkIfExists: true)
    ]

    SEQKIT_GREP_FILE ( sequence, pattern )
}

workflow test_seqkit_grep_string {

    sequence = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_gz'], checkIfExists: true)
    ]

    pattern = []

    SEQKIT_GREP_STRING ( sequence, pattern)
}
