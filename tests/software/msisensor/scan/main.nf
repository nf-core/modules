#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSOR_SCAN } from '../../../../software/msisensor/scan/main.nf' addParams( options: [:] )

workflow test_msisensor_scan {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]

    MSISENSOR_SCAN ( input )
}
