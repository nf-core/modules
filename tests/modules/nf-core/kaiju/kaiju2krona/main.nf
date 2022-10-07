#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../../../modules/nf-core/untar/main.nf'
include { KAIJU_KAIJU } from '../../../../../modules/nf-core/kaiju/kaiju/main.nf'
include { KAIJU_KAIJU2KRONA } from '../../../../../modules/nf-core/kaiju/kaiju2krona/main.nf'

workflow test_kaiju_kaiju2krona {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    db = [ [], file(params.test_data['sarscov2']['genome']['kaiju_tar_gz'], checkIfExists: true) ]

    UNTAR ( db )
    KAIJU_KAIJU ( input, UNTAR.out.untar.map{ it[1] } )
    KAIJU_KAIJU2KRONA ( KAIJU_KAIJU.out.results , UNTAR.out.untar.map{ it[1] } )
}
