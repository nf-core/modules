#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YARA } from '../../../software/yara/mapper/main.nf' addParams( options: ['args': '-e 3'] )

workflow test_yara {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]

    YARA ( input )
}
