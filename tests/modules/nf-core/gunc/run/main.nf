#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { GUNC_RUN        } from "$moduleDir/modules/nf-core/gunc/run/main.nf"
include { GUNC_DOWNLOADDB } from "$moduleDir/modules/nf-core/gunc/downloaddb/main.nf"

workflow test_gunc_run {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    GUNC_DOWNLOADDB ( 'progenomes' )
    GUNC_RUN ( input, GUNC_DOWNLOADDB.out.db )
}
