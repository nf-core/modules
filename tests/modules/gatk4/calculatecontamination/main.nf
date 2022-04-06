#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CALCULATECONTAMINATION } from '../../../../modules/gatk4/calculatecontamination/main.nf'

workflow test_gatk4_calculatecontamination_tumor_only {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test2_pileups_table'], checkIfExists: true),
              [] ]

    segmentout = false

    GATK4_CALCULATECONTAMINATION ( input, segmentout )
}

workflow test_gatk4_calculatecontamination_matched_pair {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test2_pileups_table'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_pileups_table'], checkIfExists: true) ]

    segmentout = false

    GATK4_CALCULATECONTAMINATION ( input, segmentout )
}

workflow test_gatk4_calculatecontamination_segmentation {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test2_pileups_table'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_pileups_table'], checkIfExists: true) ]

    segmentout = true

    GATK4_CALCULATECONTAMINATION ( input, segmentout )
}
