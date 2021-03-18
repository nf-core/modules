#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTP } from '../../../software/fastp/main.nf'  addParams( options: [:] )

/*
 * Test with single-end data
 */
workflow test_fastp_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ] ]
    FASTP ( input )
}

/*
 * Test with paired-end data
 */
workflow test_fastp_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true) ] ]
    FASTP ( input )
}

