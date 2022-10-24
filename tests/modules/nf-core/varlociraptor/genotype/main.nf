#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_GENOTYPE } from '../../../../../modules/nf-core/varlociraptor/genotype/main.nf'

workflow test_varlociraptor_genotype {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_GENOTYPE ( input )
}
