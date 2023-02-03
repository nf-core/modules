#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_SPLITFASTQBYLANE } from '../../../../../modules/nf-core/custom/splitfastqbylane/main.nf'

workflow test_custom_splitfastqbylane {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
	]
    ]

    CUSTOM_SPLITFASTQBYLANE ( input )
}

workflow test_custom_splitfastqbylane_single_end { 

    input = [
         [ id:'test', single_end:true ],
	 file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    CUSTOM_SPLITFASTQBYLANE ( input )
}
