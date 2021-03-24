#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YARA_INDEX } from '../../../../software/yara/index/main.nf' addParams( options: [:] )

workflow test_yara_index {

    def input = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    YARA_INDEX ( input )
}
