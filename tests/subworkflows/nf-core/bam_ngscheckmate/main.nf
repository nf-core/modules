#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_NGSCHECKMATE } from '../../../../subworkflows/nf-core/bam_ngscheckmate/main.nf'

workflow test_bam_ngscheckmate {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    BAM_NGSCHECKMATE ( input )
}
