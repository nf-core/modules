#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE_PHASE       } from '../../../../../modules/nf-core/glimpse/phase/main.nf'
include { GLIMPSE_LIGATE      } from '../../../../../modules/nf-core/glimpse/ligate/main.nf'
include { GLIMPSE_CONCORDANCE } from '../../../../../modules/nf-core/glimpse/concordance/main.nf'

workflow test_glimpse_concordance {
    
   input_vcf = [
        [ id:'input', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true),
        "chr21:16600000-16800000",
        "chr21:16650000-16750000"
    ]
    
    ref_panel = [
        [ id:'reference', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ]

    GLIMPSE_PHASE (
        input_vcf,
        ref_panel,
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true),
        []) // [meta, vcf, index, region_in, region_out], map, sample, [meta, ref, index], [meta, txt]

    all_files = GLIMPSE_PHASE.output.phased_variant
                .map { it[1] }
                .collectFile(){ item ->[ "all_files.txt", "$item" + '\n' ]}

    ligate_input=Channel.of([[ id:'input', single_end:false ]]).combine(all_files)
    GLIMPSE_LIGATE ( ligate_input )
    
    allele_freq = [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz",checkIfExists:true),
                   file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz.csi",checkIfExists:true)]
    
    truth = [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.bcf",checkIfExists:true),
                  file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.bcf.csi",checkIfExists:true)]
    
    list_inputs = Channel.of(["chr21", allele_freq[0], truth[0]])
                    .combine(GLIMPSE_LIGATE.out.merged_variants.map{it[1]}.collect().map{it[0]})
                    .collect()
    concordance_input=Channel.of([[ id:'input', single_end:false ]]).combine(list_inputs)

    GLIMPSE_CONCORDANCE ( concordance_input, [], [], []) // meta, Region, Frequencies, Truth, Estimate, minPROB, minDP, bins
}
