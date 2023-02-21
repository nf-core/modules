#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE_PHASE } from '../../../../../modules/nf-core/glimpse/phase/main.nf'

input_vcf = Channel.of([
    [ id:'input', single_end:false ], // meta map
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true),
    "chr21:16600000-16800000",
    "chr21:16650000-16750000",
    "1"
])

ref_panel = [
        [ id:'reference', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ]

samples_file = Channel
    .of('NA12878 2')
    .collectFile(name: 'sampleinfos.txt')

workflow test_glimpse_phase_nosampleinfos {
    GLIMPSE_PHASE (
        input_vcf.combine(Channel.of([[]])),
        ref_panel,
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
        ) // [meta, vcf, index, regionin, regionout, regionindex, sample_infos], map, sample, [meta, ref, index]
}

workflow test_glimpse_phase_withsampleinfos {
    GLIMPSE_PHASE (
        input_vcf.combine(samples_file),
        ref_panel,
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
        ) // [meta, vcf, index, regionin, regionout, regionindex, sample_infos], map, sample, [meta, ref, index]
}
