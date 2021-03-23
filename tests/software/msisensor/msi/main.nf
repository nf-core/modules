#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSOR_SCAN } from '../../../../software/msisensor/scan/main.nf' addParams( options: [:] )
include { MSISENSOR_MSI } from '../../../../software/msisensor/msi/main.nf' addParams( options: [:] )

workflow test_msisensor_msi {

    def scaninput  = []
    scaninput = [ [ id:'test', single_end:false ], // meta map
                  file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]
    scan = MSISENSOR_SCAN ( scaninput )

    // IMPERFECT TEST:
    // USING SARS-COV2 DATA AS NORMAL:TUMOR PAIR THIS WILL SUFFICE TO
    // TEST MODULE EXECUTION, BUT NOT FUNCTIONALITY.
    // FUNCTIONALITY HAS BEEN TESTED MANUALY USING AUTHOR-PROVIDED TEST
    // DATA (https://github.com/ding-lab/msisensor/tree/master/test)
    def input      = []

    input = Channel.from([ [ id:'test', single_end:false ], // meta map
               file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.sorted.bam", checkIfExists: true),
               file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.sorted.bam.bai", checkIfExists: true),
               file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true),
               file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) ])

    // BIT CLUMSY:
    MSISENSOR_MSI ( input.mix(scan.txt).collect() )
}
