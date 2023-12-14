#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE2_CHUNK } from '../../../../../modules/nf-core/glimpse2/chunk/main.nf'

input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true),
        "chr21"]

workflow test_glimpse2_chunk {
    ch_map_empty = Channel.of([
        [ id:'map'],
        []
    ]).collect()

    GLIMPSE2_CHUNK (input, ch_map_empty, "recursive")

}

workflow test_glimpse2_chunk_withmap {
    ch_map = Channel.of([
        [ id:'map'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
    ]).collect()

    GLIMPSE2_CHUNK (input,ch_map,"recursive")
}

