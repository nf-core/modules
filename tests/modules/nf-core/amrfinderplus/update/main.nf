#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMRFINDERPLUS_UPDATE } from '../../../../modules/amrfinderplus/update/main.nf'

workflow test_amrfinderplus_update {

    AMRFINDERPLUS_UPDATE ( )

}
