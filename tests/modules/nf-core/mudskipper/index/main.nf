#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MUDSKIPPER_INDEX } from '../../../../../modules/nf-core/mudskipper/index/main.nf'

workflow test_mudskipper_index {
    
    input = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    MUDSKIPPER_INDEX ( input )
}
