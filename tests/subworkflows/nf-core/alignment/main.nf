#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALIGNMENT } from '../../../../subworkflows/nf-core/alignment/main.nf'

workflow test_alignment {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    ALIGNMENT ( input )
}
