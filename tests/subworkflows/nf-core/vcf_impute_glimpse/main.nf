#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_IMPUTE_GLIMPSE } from '../../../../subworkflows/nf-core/vcf_impute_glimpse/main.nf'

    samples_infos = Channel.of('NA12878 2').collectFile(name: 'sampleinfos.txt')
    empty_channel = Channel.of([[]])
    region        = Channel.of(["chr21"])
    input_vcf     = Channel.of([
        [ id:'input', region: 'chr21'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true)        
    ])

    input_vcf_with_samples_infos    = input_vcf.combine(samples_infos).combine(region)
    input_vcf_without_samples_infos = input_vcf.combine(empty_channel).combine(region)

    ch_ref_panel = Channel.of([
        [ ref:'ref_panel'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true)
    ])

    ch_map = Channel.of([
        [ map:'map'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
    ])

workflow test_vcf_impute_glimpse_without_sample {
    VCF_IMPUTE_GLIMPSE ( input_vcf_without_samples_infos, ch_ref_panel, ch_map )
}

workflow test_vcf_impute_glimpse_with_sample {
    VCF_IMPUTE_GLIMPSE ( input_vcf_with_samples_infos, ch_ref_panel, ch_map )
}

workflow test_vcf_impute_glimpse_multiple {
    region        = Channel.fromList(["chr21:16600000-16790000", "chr21:16610000-16800000"])
    input_vcf     = Channel.fromList([
        [[ id:'input'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true)],
        [[ id:'input2'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true)]
    ])
    input_vcf_multiple    = input_vcf
        .combine(Channel.of([[]]).collect())
        .combine(region)
        .map{ metaI, vcf, index, sample, region ->
            [metaI + ["region": region], vcf, index, sample, region]}
        .view()
    ch_ref_panel = Channel.fromList([
        [[ ref:'ref_panel'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true)],
        [[ ref:'ref_panel2'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true)]
    ])
    VCF_IMPUTE_GLIMPSE ( input_vcf_multiple, ch_ref_panel, ch_map )
}
