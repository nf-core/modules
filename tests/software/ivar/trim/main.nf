#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_TRIM } from '../../../../software/ivar/trim/main.nf' addParams([:])

workflow test_ivar_trim {
    bed_file = file("${launchDir}/tests/data/genomics/sarscov2/bed/test-sc2-artic-v3.bed", checkIfExists: true)

    def input = []
    input = [ [ id:'test'],
                file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3.bam", checkIfExists: true) ]

    IVAR_TRIM ( input, bed_file )
}
