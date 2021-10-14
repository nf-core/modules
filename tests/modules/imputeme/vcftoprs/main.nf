#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IMPUTEME_VCFTOPRS } from '../../../../modules/imputeme/vcftoprs/main.nf' addParams( options: [:] )

workflow test_imputeme_vcftoprs {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    IMPUTEME_VCFTOPRS ( input )
}
