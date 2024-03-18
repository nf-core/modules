#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIPLE_IMPUTE_GLIMPSE2 } from '../../../../subworkflows/nf-core/multiple_impute_glimpse2/main.nf'

    ch_input_vcf = Channel.of([
        [ id:'input_vcf'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi",
            checkIfExists: true),
    ])
    sample = Channel.of('NA12878 2')
            .collectFile(name: 'sampleinfos.txt')
    ch_input_vcf_with_sample = ch_input_vcf.combine(sample)
    ch_input_vcf_without_sample = ch_input_vcf.combine(Channel.of([[]]))

    ch_ref_panel = Channel.of([
        [ id:'ref_panel'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true),
        "chr21"
    ])
    ch_map = Channel.of([
        [ id:'map'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
    ]).collect()

    ch_fasta_empty = Channel.of([
        [ id:'ref_fasta'],
        [],
        []
    ]).collect()

workflow test_multiple_impute_glimpse2_without_sample {
    MULTIPLE_IMPUTE_GLIMPSE2 ( ch_input_vcf_without_sample,
                        ch_ref_panel, ch_map, ch_fasta_empty, "recursive" )
}

workflow test_multiple_impute_glimpse2_with_sample {
    MULTIPLE_IMPUTE_GLIMPSE2 ( ch_input_vcf_with_sample,
                        ch_ref_panel, ch_map, ch_fasta_empty, "recursive" )
}

workflow test_multiple_impute_glimpse2_bam {

    input_bam = Channel.of([
        [id:'input'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.bam", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.bam.bai", checkIfExists: true),
        []
    ])

    MULTIPLE_IMPUTE_GLIMPSE2 ( input_bam, ch_ref_panel, ch_map, ch_fasta_empty, "recursive" )
}

workflow test_multiple_impute_glimpse2_cram {
    input_cram = Channel.of([
        [id:'input'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.cram", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.cram.crai", checkIfExists: true),
        []
    ])
    ch_fasta = Channel.of([
                    [id:'refHG38_chr21'],
                    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/hs38DH.chr21.fa.gz", checkIfExists: true),
                    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/hs38DH.chr21.fa.gz.fai", checkIfExists: true)
                ]).collect()
    MULTIPLE_IMPUTE_GLIMPSE2 ( input_cram, ch_ref_panel, ch_map, ch_fasta, "recursive" )
}
