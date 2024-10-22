#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_PBCSTAT } from '../../../../../modules/nf-core/purgedups/pbcstat/main.nf'
include { PURGEDUPS_CALCUTS } from '../../../../../modules/nf-core/purgedups/calcuts/main.nf'
include { PURGEDUPS_HISTPLOT } from '../../../../../modules/nf-core/purgedups/histplot/main.nf'

workflow test_purgedups_histplot {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_paf'], checkIfExists: true)
    ]

    PURGEDUPS_PBCSTAT ( input )
    PURGEDUPS_CALCUTS ( PURGEDUPS_PBCSTAT.out.stat )
    PURGEDUPS_HISTPLOT ( PURGEDUPS_PBCSTAT.out.stat.join( PURGEDUPS_CALCUTS.out.cutoff ) )
}
