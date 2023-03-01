#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE2_CHUNK } from '../../../../../modules/nf-core/glimpse2/chunk/main.nf'

workflow test_glimpse2_chunk {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
        "chr21",
        []]
    GLIMPSE2_CHUNK (input, "recursive")

}

workflow test_glimpse2_chunk_withmap {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
        "chr21",
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)]
    GLIMPSE2_CHUNK (input, "recursive")

}

