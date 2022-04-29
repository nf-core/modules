#!/usr/bin/env nextflow



include { SNAPALIGNER_INDEX } from '../../../../modules/snapaligner/index/main.nf'
include { SNAPALIGNER_SINGLE } from '../../../../modules/snapaligner/single/main.nf'

workflow test_snapaligner_single {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
    ]

    SNAPALIGNER_INDEX ( file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),[],[],[])
    SNAPALIGNER_SINGLE ( input, SNAPALIGNER_INDEX.out.index )
}
