#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_GENOTYPEGVCF } from '../../../../../modules/nf-core/parabricks/genotypegvcf/main.nf'

workflow test_parabricks_genotypegvcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PARABRICKS_GENOTYPEGVCF ( input )
}
