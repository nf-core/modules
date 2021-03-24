#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BEDTOINTERVALLIST } from '../../../../software/gatk4/bedtointervallist/main.nf' addParams( options: [:] )

workflow test_gatk4_bedtointervallist {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)] 
            ]
    dict  = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)

    GATK4_BEDTOINTERVALLIST ( input, dict )
}
