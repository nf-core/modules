#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_INDEX       } from '../../../../../modules/nf-core/bcftools/index/main.nf'
include { GLIMPSE2_PHASE       } from '../../../../../modules/nf-core/glimpse2/phase/main.nf'
include { GLIMPSE_LIGATE       } from '../../../../../modules/nf-core/glimpse/ligate/main.nf'
include { GLIMPSE2_CONCORDANCE } from '../../../../../modules/nf-core/glimpse2/concordance/main.nf'
include {
    GLIMPSE2_CONCORDANCE as GLIMPSE2_CONCORDANCE_R2_PER_SITE
                               } from '../../../../../modules/nf-core/glimpse2/concordance/main.nf'

workflow preprocessing {
    samples_infos = Channel.of('NA12878 2').collectFile(name: 'sampleinfos.txt')
    region        = Channel.of(["chr21:16600000-16800000","chr21:16650000-16750000"])
    input_vcf     = Channel.of([
        [ id:'input'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true)
    ])

    input_vcf_with_samples_infos = input_vcf.combine(samples_infos).combine(region)


    ref_panel = Channel.of([
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ])

    ch_map = Channel.of([
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true),
    ])

    GLIMPSE2_PHASE (
        input_vcf_with_samples_infos.combine( ref_panel )
                                    .combine( ch_map ),
        Channel.of([[],[],[]])
    ) // [meta, vcf, index, sample_infos, regionin, regionout, regionindex, ref, ref_index, map]

    allele_freq = Channel.of([
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz",checkIfExists:true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz.csi",checkIfExists:true)
    ])

    truth = Channel.of([
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.bcf",checkIfExists:true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.bcf.csi",checkIfExists:true)
    ])

    BCFTOOLS_INDEX ( GLIMPSE2_PHASE.output.phased_variants )

    list_inputs = GLIMPSE2_PHASE.output.phased_variants
                        .join( BCFTOOLS_INDEX.out.csi )
                        .combine( truth )
                        .combine( allele_freq )
                        .combine( Channel.of([[]]) )
                        .combine( Channel.of(["chr21"]) )

    emit:
    list_inputs
}


workflow test_glimpse2_concordance {
    list_inputs = preprocessing()
    GLIMPSE2_CONCORDANCE ( list_inputs,
                        Channel.of([[id:"params"],[],"0 0.01 0.05 0.1 0.2 0.5",[],[]]),
                        0.9999,
                        8) // [meta, estimate, estimate_index, truth, truth_index, freq, freq_index, samples, region], [meta, group, bins, ac_bins, allele_count], min-val-gl, min-val-dp
}


workflow test_list_region {
    list_inputs_list_region = preprocessing().map{
        def reg = it[8]
        def ret = it.clone()
        ret[8] = [reg, reg]
        return ( ret )
    }

    GLIMPSE2_CONCORDANCE ( list_inputs_list_region,
                        Channel.of([[id:"params"],[],"0 0.01 0.05 0.1 0.2 0.5",[],[]]),
                        0.9999,
                        8) // [meta, estimate, estimate_index, truth, truth_index, freq, freq_index, samples, region], [meta, group, bins, ac_bins, allele_count], min-val-gl, min-val-dp
}


workflow test_r2_per_site {
    list_inputs = preprocessing()
    GLIMPSE2_CONCORDANCE_R2_PER_SITE ( list_inputs,
                        Channel.of([[id:"params"],[],"0 0.01 0.05 0.1 0.2 0.5",[],[]]),
                        0.9999,
                        8) // [meta, estimate, estimate_index, truth, truth_index, freq, freq_index, samples, region], [meta, group, bins, ac_bins, allele_count], min-val-gl, min-val-dp
}
