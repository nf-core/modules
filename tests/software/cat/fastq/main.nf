#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CAT_FASTQ } from '../../../../software/cat/fastq/main.nf' addParams( options: [:] )

workflow test_cat_fastq_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
                [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                  file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.2.fastq.gz", checkIfExists: true) ]]

    CAT_FASTQ ( input )
}

workflow test_cat_fastq_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
                [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                  file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true),
                  file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.2.fastq.gz", checkIfExists: true),
                  file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.2.fastq.gz", checkIfExists: true) ]]

    CAT_FASTQ ( input )
}
