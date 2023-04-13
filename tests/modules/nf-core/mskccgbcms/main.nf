#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GBCMS } from '../../../../modules/nf-core/mskccgbcms/main.nf'

workflow test_gbcms {
    
    input = [
        [ id:'test' ], // meta map
        file(params.fasta, checkIfExists: true),
        file(params.fastafai, checkIfExists: true),
        file(params.bam, checkIfExists: true),
        file(params.bambai, checkIfExists: true),
        file("tools-test-dataset/getbasecountmultisample/chr22.maf"),
        params.sample,
        // customize test specific parameters here
        "variant_file.maf",
        [args: '--omaf']
        
    ]

    GBCMS ( input )
}

