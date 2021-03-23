#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTO_QUANT } from '../../../../software/kallisto/quant/main.nf' addParams( options: [:] )

workflow test_kallisto_quant {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]

    KALLISTO_QUANT ( input )
}
