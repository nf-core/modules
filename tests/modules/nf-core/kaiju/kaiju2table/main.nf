#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR             } from '../../../../../modules/nf-core/untar/main.nf'
include { KAIJU_KAIJU       } from '../../../../../modules/nf-core/kaiju/kaiju/main.nf'
include { KAIJU_KAIJU2TABLE } from '../../../../../modules/nf-core/kaiju/kaiju2table/main.nf'

workflow test_kaiju_kaiju_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    db    = [ [], file(params.test_data['sarscov2']['genome']['kaiju_tar_gz'], checkIfExists: true) ]
    taxon_rank = "species"

    ch_db = UNTAR ( db )
    KAIJU_KAIJU ( input, ch_db.untar.map{ it[1] } )
    KAIJU_KAIJU2TABLE ( KAIJU_KAIJU.out.results, ch_db.untar.map{ it[1] }, taxon_rank )
}
