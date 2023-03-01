#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHAPEIT5_PHASECOMMON } from '../../../../../modules/nf-core/shapeit5/phasecommon/main.nf'
include { SHAPEIT5_SWITCH } from '../../../../../modules/nf-core/shapeit5/switch/main.nf'

workflow test_shapeit5_switch {
    input_vcf = Channel.of([
    [ id:'input', single_end:false ], // meta map
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
    "chr21",
    []
    ])

    ref_panel = Channel.of([[],[],[]])
    scaffold  = Channel.of([[],[],[]])
    map       = Channel.of([[ id:'map', single_end:false ],
                            [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)]])

    SHAPEIT5_PHASECOMMON ( input_vcf, ref_panel, scaffold, map)
    
    allele_freq = [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz",checkIfExists:true),
                   file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz.csi",checkIfExists:true)]
    
    truth = [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",checkIfExists:true),
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",checkIfExists:true)]

    concordance_input=Channel.of([[ id:'input', single_end:false ]])
                            .combine(SHAPEIT5_PHASECOMMON.out.phased_variant)
                            .concat(Channel.of("chr21"))
                            .concat(truth[0])
                            .concat(allele_freq[0])
                            .concat(Channel.of([[]]))

    SHAPEIT5_SWITCH ( concordance_input) // meta, Region, Frequencies, Truth, Estimate, minPROB, minDP, bins
}
