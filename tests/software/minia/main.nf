#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIA } from '../../../../software/minia/main.nf' addParams( options: [:] )

workflow test_minia {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]
    
    

    cookiecutter.tool_name_upper ( input )
}