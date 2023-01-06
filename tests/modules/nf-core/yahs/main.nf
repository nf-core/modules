#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YAHS } from '../../../../modules/nf-core/yahs/main.nf'

workflow test_yahs {
    
    input_bed = [ [ id:'test' ], 
            file( params.test_data['yahs']['ilEupCent1_sorted_bed'], checkIfExists: true) 
            ] // meta map
    input_fa = [ file(params.test_data['yahs']['ilEupCent1_fa'], checkIfExists: true),
          file(params.test_data['yahs']['ilEupCent1_fa_fai'], checkIfExists: true)
        ]

    YAHS ( input_bed, input_fa, 'true', '', '' )
}
