#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXTMAP } from '../../../../modules/nf-core/pretextmap/main.nf'

workflow test_pretextmap_bam {
    
    input = [
        [ id:'test_bam', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    PRETEXTMAP ( input )
}

workflow test_pretextmap_cram {

    input = [
        [ id:'test_cram', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
    ]

    PRETEXTMAP ( input )
}

workflow test_pretextmap_pairs_gz {

    input = [
        [ id:'test_pairs_gz', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_pairs_gz'], checkIfExists: true)
    ]

    PRETEXTMAP ( input )
}
