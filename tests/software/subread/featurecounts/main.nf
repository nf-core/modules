#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBREAD_FEATURECOUNTS } from '../../../../software/subread/featurecounts/main.nf' addParams( options: [:] )

workflow test_subread_featurecounts_single_end {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) ]

    SUBREAD_FEATURECOUNTS ( input )
}

workflow test_subread_featurecounts_paired_end {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) ]

    SUBREAD_FEATURECOUNTS ( input )
}