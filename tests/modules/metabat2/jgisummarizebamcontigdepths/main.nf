#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../../modules/metabat2/jgisummarizebamcontigdepths/main.nf'

workflow test_metabat2_jgisummarizebamcontigdepths {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( input )
}
