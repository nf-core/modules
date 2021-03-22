#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPTITYPE_TYPE } from '../../../../software/optitype/type/main.nf' addParams( options: [:] )

workflow test_optitype_type {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]

    OPTITYPE_TYPE ( input )
}
