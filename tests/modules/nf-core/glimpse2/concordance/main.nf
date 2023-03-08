#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE2_PHASE       } from '../../../../../modules/nf-core/glimpse2/phase/main.nf'
include { GLIMPSE_LIGATE      } from '../../../../../modules/nf-core/glimpse/ligate/main.nf'
include { GLIMPSE2_CONCORDANCE } from '../../../../../modules/nf-core/glimpse2/concordance/main.nf'
include { BCFTOOLS_INDEX  } from '../../../../../modules/nf-core/bcftools/index/main.nf'

workflow test_glimpse2_concordance {
    
   input_vcf = Channel.of([
        [ id:'input', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true),
        "chr21:16600000-16800000",
        "chr21:16650000-16750000"
    ])
    
    ref_panel = Channel.of([
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ])

    ch_map = Channel.of([
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true),
    ])

    samples_file = Channel.of('NA12878 2')
                    .collectFile(name: 'sampleinfos.txt')

    GLIMPSE2_PHASE (
        input_vcf.combine(ref_panel)
                .combine(ch_map)
                .combine(samples_file),
        Channel.of([[],[],[]])
    ) // [meta, vcf, index, regionin, regionout, regionindex, ref, ref_index, map, sample_infos]
   
    allele_freq = [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz",checkIfExists:true),
                   file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz.csi",checkIfExists:true)]
    
    truth = [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.bcf",checkIfExists:true),
                  file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.bcf.csi",checkIfExists:true)]
    
    list_inputs = Channel.of(["chr21", allele_freq[0], truth[0]])
                    .combine(GLIMPSE2_PHASE.output.phased_variant.map{it[1]})
                    .combine(Channel.of([[]]))
                    .collect()
                    .view()
    concordance_input = Channel.of([[ id:'input', single_end:false ]]).combine(list_inputs)

    GLIMPSE2_CONCORDANCE ( concordance_input, Channel.of([[],[]])) // meta, Region, Frequencies, Truth, Estimate, minPROB, minDP, bins
}
