#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { AMRFINDERPLUS_UPDATE}         from '../../../../modules/amrfinderplus/update/main.nf'
include { AMRFINDERPLUS_RUN}            from '../../../../modules/amrfinderplus/run/main.nf'
include { HAMRONIZATION_AMRFINDERPLUS } from '../../../../modules/hamronization/amrfinderplus/main.nf'

workflow test_hamronization_amrfinderplus {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/amrfinderplus/test_output.tsv", checkIfExists: true)
    ]

    //AMRFINDERPLUS_UPDATE ( )
    //AMRFINDERPLUS_RUN ( input, AMRFINDERPLUS_UPDATE.out.db)
    //HAMRONIZATION_AMRFINDERPLUS ( AMRFINDERPLUS_RUN.out.report, 'tsv', AMRFINDERPLUS_RUN.out.versions, '2022-05-26.1' )
    HAMRONIZATION_AMRFINDERPLUS ( input, 'tsv', '3.10.30', '2022-05-26.1' )
}
