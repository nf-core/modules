#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_FASTQTOBAM } from '../../../../modules/fgbio/fastqtobam/main.nf' addParams( options: [:] )

workflow test_fgbio_fastqtobam {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    FGBIO_FASTQTOBAM ( input )
}
