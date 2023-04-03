#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_PHASE_SHAPEIT5 } from '../../../../subworkflows/nf-core/vcf_phase_shapeit5/main.nf'

workflow test_vcf_phase_shapeit5 {
    
    ref_panel = Channel.of([
        [ id:'input1',region:"chr21:16600115-16799989"],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true),
        [],
        "chr21:16600115-16799989"
    ])

    VCF_PHASE_SHAPEIT5 ( ref_panel, Channel.of([[],[],[]]).collect(),Channel.of([[],[],[]]).collect(),Channel.of([[],[]]).collect() )
}
