#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUTADAPT } from '../../../software/cutadapt/main.nf'  addParams( options: [ args:'-q 25' ] )

/*
 * Test with single-end data
 */
workflow test_cutadapt_se {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_single_end.fastq.gz", checkIfExists: true) ] ]

    CUTADAPT( input )
}

/*
 * Test with paired-end data
 */

workflow test_cutadapt_pe {
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/rna/test_R2.fastq.gz", checkIfExists: true) ] ]

    CUTADAPT( input )
}

