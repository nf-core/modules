#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTK_FASTK  } from '../../../modules/fastk/fastk/main.nf'
include { FASTK_HISTEX } from '../../../modules/fastk/histex/main.nf'
include { GENESCOPEFK  } from '../../../modules/genescopefk/main.nf'

workflow test_genescopefk {

    input = [
        [ id:'test' , single_end: true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    FASTK_FASTK ( input )
    FASTK_HISTEX ( FASTK_FASTK.out.hist )
    GENESCOPEFK ( FASTK_HISTEX.out.hist )
}
