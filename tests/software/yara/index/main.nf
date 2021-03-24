#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YARA_INDEX } from '../../../../software/yara/index/main.nf' addParams( options: [:] )

workflow test_yara_index {
    
    def input = file("${launchDir}/tests/data/genomics/sarscov2/bam/test_single_end.bam", checkIfExists: true)

    YARA_INDEX ( input )
}
