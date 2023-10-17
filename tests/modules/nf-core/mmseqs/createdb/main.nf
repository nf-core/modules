#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MMSEQS_CREATEDB } from '../../../../../modules/nf-core/mmseqs/createdb/main.nf'

workflow test_mmseqs_createdb {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    MMSEQS_CREATEDB ( input )
}
