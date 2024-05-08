#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { PAFTOOLS_SAM2PAF } from '../../../../../modules/nf-core/paftools/sam2paf/main.nf'

workflow test_paftools_sam2paf {
    
    input = [
        [id:'test'],
        file(params.test_data['homo_sapiens']['scramble']['bam'],
        checkIfExists: true)
    ]

    PAFTOOLS_SAM2PAF ( input )

}
