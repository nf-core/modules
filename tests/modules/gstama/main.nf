#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GSTAMA } from '../../../software/gstama/main.nf' addParams( options: [ args:"-x capped -b BAM" ] )

workflow test_gstama {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['galgal6']['bam']['chr30_refined_mapped'], checkIfExists: true)
            ]
    genome = file(params.test_data['galgal6']['fasta']['chr30'], checkIfExists: true)

    GSTAMA ( input, genome )
}
