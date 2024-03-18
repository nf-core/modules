#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HIFIASM } from '../../../../modules/nf-core/hifiasm/main.nf'

/*
 * Test with long reads only
 */
workflow test_hifiasm_hifi_only {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
    ]

    HIFIASM ( input, [], [], [], [] )
}

/*
 * Test with parental reads for phasing
 */
workflow test_hifiasm_with_parental_reads {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
    ]
    paternal_kmer_dump = file(params.test_data['homo_sapiens']['illumina']['test_yak'], checkIfExists: true)
    maternal_kmer_dump = file(params.test_data['homo_sapiens']['illumina']['test2_yak'], checkIfExists: true)

    HIFIASM ( input, paternal_kmer_dump, maternal_kmer_dump, [], [] )
}

/*
 * Test with parental reads for phasing
 */
workflow test_hifiasm_with_hic {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
    ]
    hic_reads1 = file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    hic_reads2 = file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)

    HIFIASM ( input, [], [], hic_reads1, hic_reads2 )
}

