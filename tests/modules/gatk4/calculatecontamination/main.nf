#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CALCULATECONTAMINATION } from '../../../../modules/gatk4/calculatecontamination/main.nf' addParams( options: [:] )

workflow test_gatk4_calculatecontamination_tumor_only {

    input = [ [ id:'test' ], // meta map
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test2.pileups.table", checkIfExists: true),
              [] ]

    segmentout = false

    GATK4_CALCULATECONTAMINATION ( input , segmentout )
}

workflow test_gatk4_calculatecontamination_matched_pair {

    input = [ [ id:'test' ], // meta map
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test2.pileups.table", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test.pileups.table", checkIfExists: true) ]

    segmentout = false

    GATK4_CALCULATECONTAMINATION ( input , segmentout )
}

workflow test_gatk4_calculatecontamination_segmentation {

    input = [ [ id:'test' ], // meta map
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test2.pileups.table", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test.pileups.table", checkIfExists: true) ]

    segmentout = true

    GATK4_CALCULATECONTAMINATION ( input , segmentout )
}
