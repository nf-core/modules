#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE_PHASE } from '../../../../../modules/nf-core/glimpse/phase/main.nf'

workflow test_glimpse_phase {
    input_vcf = [
        [ id:'input2', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true),
    ]
    
    ref_panel = [
        [ id:'reference', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ]

    GLIMPSE_PHASE (
        input_vcf,
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true),
        [],
        ref_panel,
        [],
        "chr21:16600000-16800000",
        "chr21:16650000-16750000") // [meta, vcf, index], map, sample, [meta, ref, index], [meta, txt], regionin, regionout
}
