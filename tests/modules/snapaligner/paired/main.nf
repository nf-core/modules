#!/usr/bin/env nextflow



include { SNAPALIGNER_INDEX } from '../../../../modules/snapaligner/index/main.nf'
include { SNAPALIGNER_PAIRED } from '../../../../modules/snapaligner/paired/main.nf'

workflow test_snapaligner_paired {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)]
    ]

    SNAPALIGNER_INDEX ( file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),[],[],[])
    SNAPALIGNER_PAIRED ( input, SNAPALIGNER_INDEX.out.index )
}
