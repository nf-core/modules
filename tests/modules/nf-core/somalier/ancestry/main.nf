#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_ANCESTRY } from '../../../../../modules/nf-core/somalier/ancestry/main.nf'

workflow test_somalier_ancestry {



    SOMALIER_ANCESTRY ( input )
}
