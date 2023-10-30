#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE2_SPLITREFERENCE } from '../../../../../modules/nf-core/glimpse2/splitreference/main.nf'

workflow test_glimpse2_splitreference_without_map {
    input = [
        [ id:'ref1000GP', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
        "chr21:16600000-16800000",
        "chr21:16600000-16800000"]
    GLIMPSE2_SPLITREFERENCE (input, [[],[]])
}
workflow test_glimpse2_splitreference_with_map {
    input = [
        [ id:'ref1000GP', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
        "chr21:16600000-16800000",
        "chr21:16600000-16800000"]
    map = [[ id:'map'],
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)]
    GLIMPSE2_SPLITREFERENCE (input, map)
}