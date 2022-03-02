#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEM_PREPAREREFERENCE    } from '../../../../modules/rsem/preparereference/main.nf'
include { RSEM_CALCULATEEXPRESSION } from '../../../../modules/rsem/calculateexpression/main.nf'

workflow test_rsem_calculateexpression {

    input = [
        [ id:'test', single_end:false, strandedness: 'forward' ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    RSEM_PREPAREREFERENCE ( fasta, gtf )
    RSEM_CALCULATEEXPRESSION( input, RSEM_PREPAREREFERENCE.out.index )
}
