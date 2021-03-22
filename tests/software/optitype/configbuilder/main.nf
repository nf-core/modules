#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPTITYPE_CONFIGBUILDER } from '../../../../software/optitype/configbuilder/main.nf' addParams( options: [:] )

workflow test_optitype_configbuilder {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]

    OPTITYPE_CONFIGBUILDER ( input )
}
