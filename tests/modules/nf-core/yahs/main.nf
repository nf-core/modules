#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YAHS } from '../../../../modules/nf-core/yahs/main.nf'

workflow test_yahs {
    
    input_bed = [ [ id:'test' ], 
            file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ] // meta map
    input_fa = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
                file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
        ]

    YAHS ( input_bed, input_fa )
}
