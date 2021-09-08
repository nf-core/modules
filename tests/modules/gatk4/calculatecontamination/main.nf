#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CALCULATECONTAMINATION } from '../../../../modules/gatk4/calculatecontamination/main.nf' addParams( options: [:] )

workflow test_gatk4_calculatecontamination {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    GATK4_CALCULATECONTAMINATION ( input )
}
