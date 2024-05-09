#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAM2LCA_ANALYZE } from '../../../../../modules/nf-core/sam2lca/analyze/main.nf'

workflow test_sam2lca_analyze {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/maxibor/sam2lca/1.1.2/tests/data/microtest.sorted.bam", checkIfExists: true),
        file("https://raw.githubusercontent.com/maxibor/sam2lca/1.1.2/tests/data/microtest.sorted.bam.bai", checkIfExists: true)
    ]

    SAM2LCA_ANALYZE ( input , [])
}
