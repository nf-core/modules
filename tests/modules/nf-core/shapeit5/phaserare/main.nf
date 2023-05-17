#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHAPEIT5_PHASECOMMON } from '../../../../../modules/nf-core/shapeit5/phasecommon/main.nf'
include { SHAPEIT5_PHASERARE   } from '../../../../../modules/nf-core/shapeit5/phaserare/main.nf'
include { BCFTOOLS_INDEX       } from '../../../../../modules/nf-core/bcftools/index/main.nf'

workflow test_shapeit5_phaserare {
    
    input_vcf = Channel.of([
    [ id:'panel', single_end:false ], // meta map
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
    [],
    "chr21",
    ])

    ref_panel = Channel.of([[],[],[]])
    scaffold  = Channel.of([[],[],[]])
    map       = Channel.of([[],[]])

    SHAPEIT5_PHASECOMMON ( input_vcf, ref_panel, scaffold, map )

    BCFTOOLS_INDEX(SHAPEIT5_PHASECOMMON.out.phased_variant.collect())

    scaffold = SHAPEIT5_PHASECOMMON.out.phased_variant
                    .join(BCFTOOLS_INDEX.out.csi)
                    .combine(Channel.of("chr21:16600000-16800000"))

    input_vcf = Channel.of([
    [ id:'input', single_end:false ], // meta map
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
    [],
    "chr21:16600000-16800000",
    ])
    SHAPEIT5_PHASERARE ( input_vcf, scaffold, map )

}
