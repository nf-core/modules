include { GLIMPSE_CHUNK                  } from '../../../modules/nf-core/glimpse/chunk/main'
include { GLIMPSE_PHASE                  } from '../../../modules/nf-core/glimpse/phase/main'
include { GLIMPSE_LIGATE                 } from '../../../modules/nf-core/glimpse/ligate/main'
include { BCFTOOLS_INDEX as INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_IMPUTE_GLIMPSE {

    take:
    ch_input      // channel (mandatory): [ meta, vcf, csi, sample, region, ref, ref_index, map ]

    main:

    ch_versions = Channel.empty()

    input_chunk = ch_input.map{
                    meta, vcf, csi, sample, region, ref, ref_index, map ->
                    [ meta, vcf, csi, region]
                }

    GLIMPSE_CHUNK ( input_chunk )
    ch_versions = ch_versions.mix( GLIMPSE_CHUNK.out.versions )

    chunk_output = GLIMPSE_CHUNK.out.chunk_chr
                                .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
                                .map { meta, it -> [meta, it["RegionIn"], it["RegionOut"]]}

    phase_input = ch_input.map{ meta, vcf, csi, sample, region, ref, ref_index, map -> [meta, vcf, csi, sample, ref, ref_index, map]}
                        .combine(chunk_output, by: 0)
                        .map{meta, vcf, csi, sample, ref, ref_index, map, regionin, regionout ->
                            [meta, vcf, csi, sample, regionin, regionout, ref, ref_index, map]}

    GLIMPSE_PHASE ( phase_input ) // [meta, vcf, index, sample_infos, regionin, regionout, ref, ref_index, map]
    ch_versions = ch_versions.mix(GLIMPSE_PHASE.out.versions )

    INDEX_PHASE ( GLIMPSE_PHASE.out.phased_variants )
    ch_versions = ch_versions.mix( INDEX_PHASE.out.versions )

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE_PHASE.out.phased_variants
        .groupTuple( by: 0 )
        .combine( INDEX_PHASE.out.csi
            .groupTuple( by: 0 ),
            by: 0
        )

    GLIMPSE_LIGATE ( ligate_input )
    ch_versions = ch_versions.mix(GLIMPSE_LIGATE.out.versions )

    INDEX_LIGATE ( GLIMPSE_LIGATE.out.merged_variants )
    ch_versions = ch_versions.mix( INDEX_LIGATE.out.versions )

    emit:
    chunk_chr              = GLIMPSE_CHUNK.out.chunk_chr           // channel: [ val(meta), txt ]
    merged_variants        = GLIMPSE_LIGATE.out.merged_variants    // channel: [ val(meta), bcf ]
    merged_variants_index  = INDEX_LIGATE.out.csi                  // channel: [ val(meta), csi ]

    versions               = ch_versions                           // channel: [ versions.yml ]
}
