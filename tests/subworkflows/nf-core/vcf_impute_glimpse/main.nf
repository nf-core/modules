#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_IMPUTE_GLIMPSE } from '../../../../subworkflows/nf-core/vcf_impute_glimpse/main.nf'

    input_vcf = Channel.of([
        [ id:'input1'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi",
            checkIfExists: true),
        "chr21"
    ])
    ref_panel = Channel.of([
        [ id:'input1'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true)
    ])
    ch_map = Channel.of([
        [ id:'input1'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
    ])

workflow test_vcf_impute_glimpse_without_sample {
    sample = Channel.of([[ id:'input1'],[]])
    VCF_IMPUTE_GLIMPSE ( input_vcf, ref_panel, ch_map, sample )
}

workflow test_vcf_impute_glimpse_with_sample {
    sample = Channel.of([[ id:'input1']])
                    .combine(Channel.of('NA12878 2')
                                .collectFile(name: 'sampleinfos.txt')
                    )
    VCF_IMPUTE_GLIMPSE ( input_vcf, ref_panel, ch_map, sample )
}
