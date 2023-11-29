#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXTSNAPSHOT as PRETEXTSNAPSHOT_ALL} from '../../../../modules/nf-core/pretextsnapshot/main.nf'
include { PRETEXTSNAPSHOT as PRETEXTSNAPSHOT_FULL} from '../../../../modules/nf-core/pretextsnapshot/main.nf'

workflow test_pretextsnapshot_all {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['galaxea_fascicularis']['hic']['pretext'], checkIfExists: true)
    ]

    PRETEXTSNAPSHOT_ALL ( input )
}

workflow test_pretextsnapshot_full {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['galaxea_fascicularis']['hic']['pretext'], checkIfExists: true)
    ]

    PRETEXTSNAPSHOT_FULL ( input )
}
