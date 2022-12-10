#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_GETSEQS } from '../../../../../modules/nf-core/purgedups/getseqs/main.nf'

workflow test_purgedups_getseqs {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    // Custom purge dups bed - params.test_data['sarscov2']['genome']['test_bed'] causes seg fault
    ch_bed = Channel.value("MT192765.1\n0\t29828\tJUNK\n").collectFile(name: 'test.bed' ).map { bed -> [ input[0], bed ]}

    PURGEDUPS_GETSEQS ( Channel.value(input).join(ch_bed) )
}
