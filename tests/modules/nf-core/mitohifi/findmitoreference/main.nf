#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOHIFI_FINDMITOREFERENCE } from '../../../../../modules/nf-core/mitohifi/findmitoreference/main.nf'

workflow test_mitohifi_findmitoreference {
    
    species = "'Phalera flavescens'"
    email = "test@email.se"
    min_length = 15659

    MITOHIFI_FINDMITOREFERENCE ( species, email, min_length )
}
