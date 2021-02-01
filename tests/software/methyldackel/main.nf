#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METHYLDACKEL_EXTRACT } from '../../../software/methyldackel/extract/main.nf' addParams( options: [:] )
include { METHYLDACKEL_MBIAS } from '../../../software/methyldackel/mbias/main.nf' addParams( options: [:] )

workflow test_methyldackel_extract {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
            file("${launchDir}/tests/data/bam/test.paired_end_methylated.sorted.bam", checkIfExists: true),
            file("${launchDir}/tests/data/bam/test.paired_end_methylated.sorted.bam.bai", checkIfExists: true),
            file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true),
            file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa.fai", checkIfExists: true) ]

    METHYLDACKEL_EXTRACT ( input )
}

workflow test_methyldackel_mbias {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
            file("${launchDir}/tests/data/bam/test.paired_end_methylated.sorted.bam", checkIfExists: true),
            file("${launchDir}/tests/data/bam/test.paired_end_methylated.sorted.bam.bai", checkIfExists: true),
            file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true),
            file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa.fai", checkIfExists: true) ]

    METHYLDACKEL_MBIAS ( input )
}
