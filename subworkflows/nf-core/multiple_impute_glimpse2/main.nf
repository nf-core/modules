include { GLIMPSE2_CHUNK                 } from '../../../modules/nf-core/glimpse2/chunk/main'
include { GLIMPSE2_SPLITREFERENCE        } from '../../../modules/nf-core/glimpse2/splitreference/main'
include { GLIMPSE2_PHASE                 } from '../../../modules/nf-core/glimpse2/phase/main'
include { GLIMPSE2_LIGATE                } from '../../../modules/nf-core/glimpse2/ligate/main'
include { BCFTOOLS_INDEX as INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow MULTIPLE_IMPUTE_GLIMPSE2 {

    take:
    ch_input    // channel (mandatory): [ meta, vcf, csi, infos ]
    ch_ref      // channel (mandatory): [ meta, vcf, csi, region ]
    ch_map      // channel  (optional): [ meta, map ]
    ch_fasta    // channel  (optional): [ meta, fasta, index ]
    chunk_model // string: model used to chunk the reference panel

    main:

    ch_versions = Channel.empty()

    // Chunk reference panel
    GLIMPSE2_CHUNK ( ch_ref, ch_map, chunk_model )
    ch_versions = ch_versions.mix( GLIMPSE2_CHUNK.out.versions.first() )

    chunk_output = GLIMPSE2_CHUNK.out.chunk_chr
                                .splitCsv(header: ['ID', 'Chr', 'RegionBuf', 'RegionCnk', 'WindowCm',
                                            'WindowMb', 'NbTotVariants', 'NbComVariants'],
                                        sep: "\t", skip: 0)
                                .map { meta, it -> [meta, it["RegionBuf"], it["RegionCnk"]]}

    // Split reference panel in bin files
    split_input = ch_ref.map{ meta, ref, index, region -> [meta, ref, index]}
                        .combine(chunk_output, by: 0)

    GLIMPSE2_SPLITREFERENCE( split_input, ch_map )
    ch_versions = ch_versions.mix( GLIMPSE2_SPLITREFERENCE.out.versions.first() )

    phase_input = ch_input.combine( GLIMPSE2_SPLITREFERENCE.out.bin_ref )
                        .map{ input_meta, input_file, input_index, input_infos,
                            panel_meta, panel_bin ->
                            [input_meta, input_file, input_index, input_infos,
                            [], [], panel_bin, [], []]
                    }/* Remove unnecessary meta maps
                        add null index as we use a bin file,
                        add null value for input and output region as we use a bin file */

    // Phase input files for each reference bin files + indexing
    GLIMPSE2_PHASE ( phase_input, ch_fasta ) // [meta, vcf, index, sample_infos, regionin, regionout, regionindex, ref, ref_index, map], [ meta, fasta, index ]
    ch_versions = ch_versions.mix( GLIMPSE2_PHASE.out.versions.first() )

    INDEX_PHASE ( GLIMPSE2_PHASE.out.phased_variant )
    ch_versions = ch_versions.mix( INDEX_PHASE.out.versions.first() )

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE2_PHASE.out.phased_variant
                                    .groupTuple()
                                    .combine( INDEX_PHASE.out.csi
                                            .groupTuple()
                                            .collect(), by: 0 )

    GLIMPSE2_LIGATE ( ligate_input )
    ch_versions = ch_versions.mix( GLIMPSE2_LIGATE.out.versions.first() )

    INDEX_LIGATE ( GLIMPSE2_LIGATE.out.merged_variants )
    ch_versions = ch_versions.mix( INDEX_LIGATE.out.versions.first() )

    emit:
    chunk_chr              = GLIMPSE2_CHUNK.out.chunk_chr           // channel: [ val(meta), txt ]
    merged_variants        = GLIMPSE2_LIGATE.out.merged_variants    // channel: [ val(meta), bcf ]
    merged_variants_index  = INDEX_LIGATE.out.csi                   // channel: [ val(meta), csi ]

    versions               = ch_versions                            // channel: [ versions.yml ]
}
