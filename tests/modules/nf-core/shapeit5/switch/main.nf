#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHAPEIT5_PHASECOMMON } from '../../../../../modules/nf-core/shapeit5/phasecommon/main.nf'
include { SHAPEIT5_SWITCH } from '../../../../../modules/nf-core/shapeit5/switch/main.nf'
include { BCFTOOLS_INDEX  } from '../../../../../modules/nf-core/bcftools/index/main.nf'

workflow test_shapeit5_switch {
    input_vcf = Channel.of([
    [ id:'input', single_end:false ], // meta map
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
    [],
    "chr21"
    ])

    ref_panel = Channel.of([[],[],[]])
    scaffold  = Channel.of([[],[],[]])
    map       = Channel.of([[ id:'map'],
                            [file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)]])

    SHAPEIT5_PHASECOMMON ( input_vcf, ref_panel, scaffold, map)
    
    allele_freq = Channel.of([[ id:'freq_file'],
                                file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz",checkIfExists:true),
                                file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.sites.vcf.gz.csi",checkIfExists:true)])
    
    truth = Channel.of([[ id:'truth_panel'],
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",checkIfExists:true),
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",checkIfExists:true)])

    BCFTOOLS_INDEX ( SHAPEIT5_PHASECOMMON.output.phased_variant )

    concordance_input=SHAPEIT5_PHASECOMMON.output.phased_variant.groupTuple()
                        .join(BCFTOOLS_INDEX.out.csi.groupTuple())
                        .combine(Channel.of("chr21"))
                        .combine(Channel.of([[]]))

    SHAPEIT5_SWITCH ( concordance_input, truth, allele_freq) // meta, Region, Frequencies, Truth, Estimate, minPROB, minDP, bins
}
