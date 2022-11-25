#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARCASHLA_EXTRACT }                        from '../../../../../modules/nf-core/arcashla/extract/main.nf'

workflow test_arcashla_extract {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    ARCASHLA_EXTRACT    ( input )

}

workflow test_arcashla_extract_single_end {

    input = [
        [ id:'test_single_end', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test3_single_end_markduplicates_sorted_bam'], checkIfExists: true)
    ]

    ARCASHLA_EXTRACT    ( input )

}
