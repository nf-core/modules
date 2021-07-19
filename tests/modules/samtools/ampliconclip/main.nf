#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_AMPLICONCLIP } from '../../../../modules/samtools/ampliconclip/main.nf' addParams( options: [:] )

workflow test_samtools_ampliconclip {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    SAMTOOLS_AMPLICONCLIP ( input )
}
