#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GSTAMA_COLLAPSE } from '../../../../software/gstama/collapse/main.nf' addParams( options: [ args:"-x capped -b BAM" ] )

workflow test_gstama_collapse {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['galgal6']['bam']['chr30_refined_mapped'], checkIfExists: true)
            ]
    genome = file(params.test_data['galgal6']['fasta']['chr30'], checkIfExists: true)

    GSTAMA_COLLAPSE ( input, genome )
}
