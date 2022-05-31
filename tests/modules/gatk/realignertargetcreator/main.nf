#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_REALIGNERTARGETCREATOR } from '../../../../modules/gatk/realignertargetcreator/main.nf'

workflow test_gatk_realignertargetcreator {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GATK_REALIGNERTARGETCREATOR ( input )
}
