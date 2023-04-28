#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE2_PHASE } from '../../../../../modules/nf-core/glimpse2/phase/main.nf'


input_vcf = Channel.of([
    [ id:'input' ], // meta map
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true),
    [],
    "chr21:16600000-16800000",
    "chr21:16650000-16750000"
    ])
input_bam = Channel.of([
    [id:'input'],
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.bam", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.bam.bai", checkIfExists: true),
    [],
    "chr21:16600000-16800000",
    "chr21:16650000-16750000",
])
input_cram = Channel.of([
    [id:'input'],
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.cram", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.cram.crai", checkIfExists: true),
    [],
    "chr21:16600000-16800000",
    "chr21:16650000-16750000",
])
ref_panel = Channel.of([
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ])

map_file = Channel.of([
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
    ])

reference_genome = Channel.of([
                    [id:'refHG38_chr21'],
                    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/hs38DH.chr21.fa.gz", checkIfExists: true),
                    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/hs38DH.chr21.fa.gz.fai", checkIfExists: true)
                    ])

workflow test_glimpse2_phase_vcf {
    GLIMPSE2_PHASE (
        input_vcf.combine(ref_panel)
                .combine(map_file),
        Channel.of([[],[],[]])
        ) // [meta, vcf, index, regionin, regionout, regionindex, sample_infos], map, sample, [meta, ref, index]
}

workflow test_glimpse2_phase_bam {
    GLIMPSE2_PHASE (
        input_bam.combine(ref_panel)
                .combine(map_file),
        Channel.of([[],[],[]])
        ) // [meta, vcf, index, regionin, regionout, regionindex, sample_infos], map, sample, [meta, ref, index]
}

workflow test_glimpse2_phase_cram {
    GLIMPSE2_PHASE (
        input_cram.combine(ref_panel)
                .combine(map_file),
        reference_genome
        ) // [meta, vcf, index, regionin, regionout, regionindex, sample_infos], map, sample, [meta, ref, index]
}
