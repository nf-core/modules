#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOHIFI_FINDMITOREFERENCE } from '../../../../../modules/nf-core/mitohifi/findmitoreference/main.nf'

workflow test_mitohifi_findmitoreference {

    input = [
        [ id:'test' ],       // meta map
        "Phalera flavescens" // species
    ]

    MITOHIFI_FINDMITOREFERENCE ( input )
}
