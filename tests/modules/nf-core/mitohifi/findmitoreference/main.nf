#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOHIFI_FINDMITOREFERENCE } from '../../../../../modules/nf-core/mitohifi/findmitoreference/main.nf'

workflow test_mitohifi_findmitoreference {
    
    species = "Homo sapiens"
    email = "verena.kutschera@scilifelab.se"
    min_length = 16000

    MITOHIFI_FINDMITOREFERENCE ( species, email, min_length )
}
