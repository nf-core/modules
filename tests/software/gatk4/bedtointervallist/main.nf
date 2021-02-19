#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BEDTOINTERVALLIST } from '../../../../software/gatk4/bedtointervallist/main.nf' addParams( options: [:] )

workflow test_gatk4_bedtointervallist {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/bed/sarscov2.bed", checkIfExists: true)] ]

    sd = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.dict", checkIfExists: true)

    GATK4_BEDTOINTERVALLIST ( input, sd )
}
