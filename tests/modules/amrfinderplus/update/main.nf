#!/usr/bin/env nextflow



include { AMRFINDERPLUS_UPDATE } from '../../../../modules/amrfinderplus/update/main.nf'

workflow test_amrfinderplus_update {

    AMRFINDERPLUS_UPDATE ( )

}
