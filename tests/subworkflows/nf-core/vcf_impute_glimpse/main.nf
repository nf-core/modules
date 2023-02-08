#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_IMPUTE_GLIMPSE } from '../../../../subworkflows/nf-core/vcf_impute_glimpse/main.nf'

workflow test_vcf_impute_glimpse {
    
    input_vcf = Channel.of([
        [ id:'input', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi",
            checkIfExists: true),
        "chr21"
    ])
    ref_panel = [
        [ id:'reference', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true)
    ]
    map = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)

    VCF_IMPUTE_GLIMPSE ( input_vcf, ref_panel, map, [] )
}
