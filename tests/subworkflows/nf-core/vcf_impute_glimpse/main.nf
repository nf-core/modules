#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_IMPUTE_GLIMPSE } from '../../../../subworkflows/nf-core/vcf_impute_glimpse/main.nf'

samples_infos = Channel.of('NA12878 2').collectFile(name: 'sampleinfos.txt')

ch_panel = Channel.fromList([
    [[ ref:'ref_panel'],
    file("https://github.com/nf-core/test-datasets/raw/imputation/data/panel/both/1000GP.chr21_22.noNA12878.s.bcf",
        checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/imputation/data/panel/both/1000GP.chr21_22.noNA12878.s.bcf.csi",
        checkIfExists: true)],
    [[ ref:'ref_panel2'],
    file("https://github.com/nf-core/test-datasets/raw/imputation/data/panel/both/1000GP.chr21_22.noNA12878.s.bcf",
        checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/imputation/data/panel/both/1000GP.chr21_22.noNA12878.s.bcf.csi",
        checkIfExists: true)]
])

workflow test_vcf_impute_glimpse_multiple_chr {

    region = Channel.fromList([
        [[chr: "chr21", region: "chr21:16600000-16800000"], "chr21:16600000-16800000"],
        [[chr: "chr22", region: "chr22:16600000-16800000"], "chr22:16600000-16800000"]
    ])

    input_vcf = Channel.fromList([
        [[ id:'input'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/imputation/data/NA12878/both/NA12878.chr21_22.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/imputation/data/NA12878/both/NA12878.chr21_22.s.1x.vcf.gz.csi", checkIfExists: true),
        ],
        [[ id:'input2'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/imputation/data/NA12878/both/NA12878.chr21_22.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/imputation/data/NA12878/both/NA12878.chr21_22.s.1x.vcf.gz.csi", checkIfExists: true),
        ]
    ])

    input_vcf_multiple = input_vcf
        .combine( samples_infos )
        .combine( region )
        .map{ metaI, vcf, index, sample, metaCR, region ->
            [metaI + metaCR, vcf, index, sample, region ]
        }

    ch_map = Channel.fromList([
        [[ chr: "chr21"],
            file("https://github.com/nf-core/test-datasets/raw/imputation/data/genetic_maps.b38/chr21.b38.gmap.gz", checkIfExists: true)
        ],
        [[ chr: "chr22"],
            file("https://github.com/nf-core/test-datasets/raw/imputation/data/genetic_maps.b38/chr22.b38.gmap.gz", checkIfExists: true)
        ]
    ])

    // Combine input and map depending on chromosome name
    ch_input_map = input_vcf_multiple
        .map{ metaIRC, vcf, index, sample, region ->
            [metaIRC.subMap(["chr"]), metaIRC, vcf, index, sample, region]
        }
        .combine(ch_map, by: 0)
        .map{ metaC, metaIRC, vcf, index, sample, region, map ->
            [metaIRC, vcf, index, sample, region, map] }

    // Combine input and map to reference panel could also be done by chromosome
    ch_input = ch_input_map
        .combine(ch_panel)
        .map{ metaIRC, vcf, index, sample, region, map, metaP, ref, ref_index ->
            [ metaIRC + metaP, vcf, index, sample, region, ref, ref_index, map ]
        }
    VCF_IMPUTE_GLIMPSE ( ch_input )
}

workflow test_vcf_impute_glimpse_multiple_region {
    region        = Channel.fromList(["chr21:16600000-16790000", "chr21:16610000-16800000"])
    input_vcf     = Channel.fromList([
        [[ id:'input'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true)],
        [[ id:'input2'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true)]
    ])
    input_vcf_multiple    = input_vcf
        .combine(Channel.of([[]]).collect())
        .combine( region )
        .map{ metaI, vcf, index, sample, region ->
            [metaI + ["region": region], vcf, index, sample, region]
        }
    
    ch_map = Channel.of([
        [ map:'map'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
    ])

    ch_input = input_vcf_multiple
        .combine(ch_panel)
        .combine(ch_map)
        .map { metaIR, vcf, index, sample, region, metaP, ref, ref_index, metaM, map ->
            [ metaIR + metaP + metaM, vcf, index, sample, region, ref, ref_index, map ]
        }
    VCF_IMPUTE_GLIMPSE ( ch_input )
}
