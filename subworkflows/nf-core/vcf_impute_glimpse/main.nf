include { GLIMPSE_CHUNK                           } from '../../../modules/nf-core/glimpse/chunk/main'
include { GLIMPSE_PHASE                           } from '../../../modules/nf-core/glimpse/phase/main'
include { GLIMPSE_LIGATE                          } from '../../../modules/nf-core/glimpse/ligate/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_IMPUTE_GLIMPSE {
    take:
    ch_vcf    // channel (mandatory): [ meta, vcf, csi, infos ]
    ch_ref    // channel (mandatory): [ meta, vcf, csi, region ]
    ch_chunks // channel (optional) : [ meta, regionin, regionout ]
    ch_map    // channel (optional) : [ meta, map ]
    chunk     // val (optional)     : boolean to activate/deactivate chunking step

    main:

    ch_versions = channel.empty()

    if (chunk == true) {
        // Error if pre-defined chunks are provided when chunking is activated
        ch_chunks
            .filter { _meta, regionin, regionout -> regionin.size() == 0 && regionout.size() == 0 }
            .ifEmpty { error("ERROR: Cannot provide pre-defined chunks (regionin) when chunk=true. Please either set chunk=false to use provided chunks, or remove input chunks to enable automatic chunking.") }

        GLIMPSE_CHUNK(ch_ref)
        ch_versions = ch_versions.mix(GLIMPSE_CHUNK.out.versions.first())

        ch_chunks = GLIMPSE_CHUNK.out.chunk_chr
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { meta, it -> [meta, it["RegionIn"], it["RegionOut"]] }
    }

    ch_chunks
        .filter { _meta, regionin, regionout -> regionin.size() > 0 && regionout.size() > 0 }
        .ifEmpty { error("ERROR: ch_chunks channel is empty. Please provide a valid channel or set chunk parameter to true.") }

    ch_chunks_panel_map = ch_chunks
        .combine(ch_ref, by: 0)
        .combine(ch_map, by: 0)

    ch_chunks_panel_map.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    phase_input = ch_vcf
        .combine(ch_chunks_panel_map)
        .map{ metaI, vcf, csi, sample, metaPC, regionin, regionout, ref, ref_index, _region, map -> [
            metaI + metaPC + ["regionout": regionout],
            vcf, csi, sample, // target input
            regionin, regionout, // chunks
            ref, ref_index, map // reference panel
        ]}

    GLIMPSE_PHASE(phase_input)
    ch_versions = ch_versions.mix(GLIMPSE_PHASE.out.versions.first())

    BCFTOOLS_INDEX_PHASE(GLIMPSE_PHASE.out.phased_variants)

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE_PHASE.out.phased_variants
        .join(
            BCFTOOLS_INDEX_PHASE.out.tbi.mix(BCFTOOLS_INDEX_PHASE.out.csi),
            failOnMismatch: true,
            failOnDuplicate: true,
        )
        .map { meta, vcf, index ->
            def keysToKeep = meta.keySet() - 'regionout'
            [meta.subMap(keysToKeep), vcf, index]
        }
        .groupTuple()

    GLIMPSE_LIGATE(ligate_input)
    ch_versions = ch_versions.mix(GLIMPSE_LIGATE.out.versions.first())

    BCFTOOLS_INDEX_LIGATE(GLIMPSE_LIGATE.out.merged_variants)

    // Join imputed and index files
    ch_vcf_index = GLIMPSE_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX_LIGATE.out.tbi.mix(BCFTOOLS_INDEX_LIGATE.out.csi),
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    chunks    = ch_chunks    // channel: [ val(meta), regionin, regionout ]
    vcf_index = ch_vcf_index // channel: [ val(meta), vcf, index ]
    versions  = ch_versions  // channel: [ versions.yml ]
}
