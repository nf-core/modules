#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_PBCSTAT   } from '../../../../../modules/nf-core/purgedups/pbcstat/main.nf'
include { PURGEDUPS_CALCUTS   } from '../../../../../modules/nf-core/purgedups/calcuts/main.nf'
include { PURGEDUPS_PURGEDUPS } from '../../../../../modules/nf-core/purgedups/purgedups/main.nf'

workflow test_purgedups_purgedups {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_paf'], checkIfExists: true)
    ]

    PURGEDUPS_PBCSTAT ( input )
    PURGEDUPS_CALCUTS ( PURGEDUPS_PBCSTAT.out.stat )
    PURGEDUPS_PURGEDUPS (
        PURGEDUPS_PBCSTAT.out.basecov
            .join( PURGEDUPS_CALCUTS.out.cutoff )
            .join( Channel.value( input ) )
    )
}
