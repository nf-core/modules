#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CALCULATECONTAMINATION } from '../../../../modules/gatk4/calculatecontamination/main.nf' addParams( options: [:] )

workflow test_gatk4_calculatecontamination_tumor_only {

    input = [ [ id:'test' ], // meta map
              file("/home/AD/gmackenz/test_data/pileupsummaries/test2.pileups.table", checkIfExists: true),
              [] ]

    GATK4_CALCULATECONTAMINATION ( input )
}

workflow test_gatk4_calculatecontamination_matched_pair {

    input = [ [ id:'test' ], // meta map
              file("/home/AD/gmackenz/test_data/pileupsummaries/test2.pileups.table", checkIfExists: true),
              file("/home/AD/gmackenz/test_data/pileupsummaries/test.pileups.table", checkIfExists: true) ]

    GATK4_CALCULATECONTAMINATION ( input )
}
