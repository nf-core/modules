#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GETORGANELLE_CONFIG } from '../../../../../modules/nf-core/getorganelle/config/main.nf'

workflow test_getorganelle_config {

    input = Channel.of( 'all', 'embplant_pt', 'embplant_mt', 'embplant_nr', 'fungus_mt', 'fungus_nr', 'animal_mt', 'other_pt' )

    GETORGANELLE_CONFIG ( input )
}
