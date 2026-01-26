include { GLIMPSE2_CHUNK                          } from '../../../modules/nf-core/glimpse2/chunk/main'
include { GLIMPSE2_SPLITREFERENCE                 } from '../../../modules/nf-core/glimpse2/splitreference/main'
include { GLIMPSE2_PHASE                          } from '../../../modules/nf-core/glimpse2/phase/main'
include { GLIMPSE2_LIGATE                         } from '../../../modules/nf-core/glimpse2/ligate/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow BAM_VCF_IMPUTE_GLIMPSE2 {
    take:
    ch_input // channel (mandatory): [ meta, vcf, csi, list, infos ]
    ch_ref // channel (mandatory): [ meta, vcf, csi, region ]
    ch_chunks // channel (optional) : [ meta, regionin, regionout ]
    ch_map // channel (optional) : [ meta, map ]
    ch_fasta // channel (optional) : [ meta, fasta, index ]
    chunk // val (optional)     : boolean to activate/deactivate chunking step
    chunk_model // val (optional)     : model file for chunking
    splitreference // val (optional)     : boolean to activate/deactivate split reference step

    main:

    ch_versions = channel.empty()

    if (chunk == true) {
        // Error if pre-defined chunks are provided when chunking is activated
        ch_chunks
            .filter { _meta, regionin, regionout -> regionin.size() == 0 && regionout.size() == 0 }
            .ifEmpty {
                error("ERROR: Cannot provide pre-defined chunks (regionin) when chunk=true. Please either set chunk=false to use provided chunks, or remove input chunks to enable automatic chunking.")
            }

        // Chunk reference panel
        ch_ref_map = ch_ref.combine(ch_map, by: 0)
        GLIMPSE2_CHUNK(ch_ref_map, chunk_model)
        ch_versions = ch_versions.mix(GLIMPSE2_CHUNK.out.versions.first())

        ch_chunks = GLIMPSE2_CHUNK.out.chunk_chr
            .splitCsv(
                header: [
                    'ID',
                    'Chr',
                    'RegionBuf',
                    'RegionCnk',
                    'WindowCm',
                    'WindowMb',
                    'NbTotVariants',
                    'NbComVariants',
                ],
                sep: "\t",
                skip: 0,
            )
            .map { meta, it -> [meta, it["RegionBuf"], it["RegionCnk"]] }
    }

    ch_chunks
        .filter { _meta, regionin, regionout -> regionin.size() > 0 && regionout.size() > 0 }
        .ifEmpty { error("ERROR: ch_chunks channel is empty. Please provide a valid channel or set chunk parameter to true.") }

    if (splitreference == true) {
        // Split reference panel in bin files
        split_input = ch_ref
            .combine(ch_chunks, by: 0)
            .combine(ch_map, by: 0)
            .map { meta, ref, index, _region, regionin, regionout, gmap ->
                [
                    meta + ["regionin": regionin, "regionout": regionout],
                    ref,
                    index,
                    regionin,
                    regionout,
                    gmap,
                ]
            }

        GLIMPSE2_SPLITREFERENCE(split_input)
        ch_versions = ch_versions.mix(GLIMPSE2_SPLITREFERENCE.out.versions.first())

        // Everything is provided by the bin file so no additional file
        ch_chunks_panel_map = GLIMPSE2_SPLITREFERENCE.out.bin_ref.map { meta, bin_ref -> [meta, [], [], bin_ref, [], []] }
    }
    else {
        ch_chunks_panel_map = ch_chunks
            .combine(ch_ref, by: 0)
            .combine(ch_map, by: 0)
            .map { meta, regionin, regionout, ref, ref_index, _region, gmap ->
                [
                    meta + ["regionin": regionin, "regionout": regionout],
                    regionin,
                    regionout,
                    ref,
                    ref_index,
                    gmap,
                ]
            }
    }

    ch_chunks_panel_map.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    ch_phase_input = ch_input
        .combine(ch_chunks_panel_map)
        .map { metaI, input, index, list, infos, metaCPM, regionin, regionout, panel, panel_index, gmap ->
            [
                metaI + metaCPM,
                input,
                index,
                list,
                infos,
                regionin,
                regionout,
                panel,
                panel_index,
                gmap,
            ]
        }

    // Impute with Glimpse2
    GLIMPSE2_PHASE(ch_phase_input, ch_fasta)
    ch_versions = ch_versions.mix(GLIMPSE2_PHASE.out.versions.first())

    // Index phased file
    BCFTOOLS_INDEX_PHASE(GLIMPSE2_PHASE.out.phased_variants)

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE2_PHASE.out.phased_variants
        .join(
            BCFTOOLS_INDEX_PHASE.out.tbi.mix(BCFTOOLS_INDEX_PHASE.out.csi),
            failOnMismatch: true,
            failOnDuplicate: true,
        )
        .map { meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionin', 'regionout']
            [meta.subMap(keysToKeep), vcf, index]
        }
        .groupTuple()

    GLIMPSE2_LIGATE(ligate_input)
    ch_versions = ch_versions.mix(GLIMPSE2_LIGATE.out.versions.first())

    BCFTOOLS_INDEX_LIGATE(GLIMPSE2_LIGATE.out.merged_variants)

    // Join imputed and index files
    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX_LIGATE.out.tbi.mix(BCFTOOLS_INDEX_LIGATE.out.csi),
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    chunks    = ch_chunks // channel: [ val(meta), regionin, regionout ]
    vcf_index = ch_vcf_index // channel: [ val(meta), vcf, csi ]
    versions  = ch_versions // channel: [ versions.yml ]
}
