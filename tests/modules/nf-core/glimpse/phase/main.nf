#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE_PHASE } from '../../../../../modules/nf-core/glimpse/phase/main.nf'

    input_vcf = Channel.of([
        [ id:'input', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true),
        "chr21:16600000-16800000",
        "chr21:16650000-16750000",
    ])

    ref_panel = Channel.of([
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ])

    ch_map = Channel.of([
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true),
    ])

    samples_file = Channel.of('NA12878 2').collectFile(name: 'sampleinfos.txt')

workflow test_glimpse_phase_nosampleinfos {
    GLIMPSE_PHASE (
        input_vcf.combine(ref_panel)
                .combine(ch_map)
                .combine(Channel.of([[]]))
    ) // [meta, vcf, index, regionin, regionout, regionindex, ref, ref_index, map, sample_infos]
}

workflow test_glimpse_phase_withsampleinfos {
    GLIMPSE_PHASE (
        input_vcf.combine(ref_panel)
                .combine(ch_map)
                .combine(samples_file)
    ) // [meta, vcf, index, regionin, regionout, regionindex, ref, ref_index, map, sample_infos]
}
