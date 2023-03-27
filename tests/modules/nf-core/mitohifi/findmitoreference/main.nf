#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOHIFI_FINDMITOREFERENCE } from '../../../../../modules/nf-core/mitohifi/findmitoreference/main.nf'

workflow test_mitohifi_findmitoreference {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    MITOHIFI_FINDMITOREFERENCE ( input )
}
