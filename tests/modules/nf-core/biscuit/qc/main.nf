#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_INDEX  } from '../../../../modules/biscuit/index/main.nf'
include { BISCUIT_QC } from '../../../../modules/biscuit/qc/main.nf'

workflow test_biscuit_qc {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX( fasta )
    BISCUIT_QC ( input, BISCUIT_INDEX.out.index )
}
