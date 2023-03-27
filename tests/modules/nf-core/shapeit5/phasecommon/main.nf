#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHAPEIT5_PHASECOMMON } from '../../../../../modules/nf-core/shapeit5/phasecommon/main.nf'

workflow test_shapeit5_phasecommon_without_map {
    
    input_vcf = Channel.of([
    [ id:'input', single_end:false ], // meta map
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
    "chr21",
    []
    ])

    ref_panel = Channel.of([[],[],[]])
    scaffold  = Channel.of([[],[],[]])
    map       = Channel.of([[],[]])

    SHAPEIT5_PHASECOMMON ( input_vcf, ref_panel, scaffold, map )
}

workflow test_shapeit5_phasecommon_with_map {
    
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
}
