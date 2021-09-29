#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GTDBTK_CLASSIFYWF } from '../../../../modules/gtdbtk/classifywf/main.nf' addParams( options: [:] )

workflow test_gtdbtk_classifywf {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    GTDBTK_CLASSIFYWF ( input )
}
