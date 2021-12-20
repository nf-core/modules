#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBREAD_FEATURECOUNTS } from '../../../../modules/subread/featurecounts/main.nf'

workflow test_subread_featurecounts_forward {
    
    def input = []
    input = [ [ id:'test', single_end:true, strandedness:'forward' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) ]

    SUBREAD_FEATURECOUNTS ( input )
}

workflow test_subread_featurecounts_reverse {
    
    def input = []
    input = [ [ id:'test', single_end:true, strandedness:'reverse' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) ]

    SUBREAD_FEATURECOUNTS ( input )
}

workflow test_subread_featurecounts_unstranded {
    
    def input = []
    input = [ [ id:'test', single_end:true, strandedness:'unstranded' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) ]

    SUBREAD_FEATURECOUNTS ( input )
}