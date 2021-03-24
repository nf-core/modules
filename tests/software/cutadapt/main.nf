#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUTADAPT } from '../../../software/cutadapt/main.nf'  addParams( options: [ args:'-q 25' ] )

/*
 * Test with single-end data
 */
workflow test_cutadapt_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true) ] 
            ]

    CUTADAPT ( input )
}

/*
 * Test with paired-end data
 */

workflow test_cutadapt_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ] 
            ]

    CUTADAPT ( input )
}

