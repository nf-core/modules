#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DELLY_CALL } from '../../../../software/delly/call/main.nf' addParams( options: [:] )

workflow test_delly_call {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true) ]
    ]

    bai = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true)
    fasta = file("${launchDir}/temp_test_data/genome.fasta", checkIfExists: true)
    fai = file("${launchDir}/temp_test_data/genome.fasta.fai", checkIfExists: true)

    DELLY_CALL ( input, bai, fasta, fai )
}
