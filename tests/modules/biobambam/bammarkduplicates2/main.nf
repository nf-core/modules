#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BIOBAMBAM_BAMMARKDUPLICATES2 } from '../../../../modules/biobambam/bammarkduplicates2/main.nf'

workflow test_biobambam_bammarkduplicates2 {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    BIOBAMBAM_BAMMARKDUPLICATES2 ( input )
}
