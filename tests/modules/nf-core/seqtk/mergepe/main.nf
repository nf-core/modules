#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { SEQTK_MERGEPE } from "$moduleDir/modules/nf-core/seqtk/mergepe/main.nf"

//
// Test with single-end data
//

workflow test_seqtk_mergepe_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]

    SEQTK_MERGEPE ( input )
}

//
// Test with paired-end data
//

workflow test_seqtk_mergepe_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    SEQTK_MERGEPE ( input )
}
