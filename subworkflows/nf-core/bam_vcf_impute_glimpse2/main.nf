include { GLIMPSE2_CHUNK                     } from '../../../modules/nf-core/glimpse2/chunk/main'
include { GLIMPSE2_SPLITREFERENCE            } from '../../../modules/nf-core/glimpse2/splitreference/main'
include { GLIMPSE2_PHASE                     } from '../../../modules/nf-core/glimpse2/phase/main'
include { GLIMPSE2_LIGATE                    } from '../../../modules/nf-core/glimpse2/ligate/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow BAM_VCF_IMPUTE_GLIMPSE2 {

    take:
    ch_input       // channel (mandatory): [ meta, vcf, csi, infos ]
    ch_ref         // channel (mandatory): [ meta, vcf, csi ]
    ch_chunks      // channel (optional) : [ meta, regionin, regionout ]
    ch_map         // channel (optional) : [ meta, map ]
    ch_fasta       // channel (optional) : [ meta, fasta, index ]
    chunk          // val (optional)     : boolean to activate/deactivate chunking step
    chunk_model    // val (optional)     : model file for chunking
    splitreference // val (optional)     : boolean to activate/deactivate split reference step

    main:

    ch_versions = channel.empty()

    if ( chunk == true ){
        // Chunk reference panel
        ch_ref_map = ch_ref.combine(ch_map, by: 0)
        GLIMPSE2_CHUNK ( ch_ref_map, chunk_model )
        ch_versions = ch_versions.mix( GLIMPSE2_CHUNK.out.versions.first() )

        ch_chunks = GLIMPSE2_CHUNK.out.chunk_chr
            .splitCsv(header: [
                'ID', 'Chr', 'RegionBuf', 'RegionCnk', 'WindowCm',
                'WindowMb', 'NbTotVariants', 'NbComVariants'
            ], sep: "\t", skip: 0)
            .map { meta, it -> [meta, it["RegionBuf"], it["RegionCnk"]]}
    }

    ch_chunks.ifEmpty{
        error "BAM_VCF_IMPUTE_GLIMPSE2: ch_chunks channel is empty. Please provide a valid channel or set chunk parameter to true."
    }

    if ( splitreference == true ) {
        // Split reference panel in bin files
        split_input = ch_ref
            .map{ meta, ref, index, _region -> [meta, ref, index] }
            .combine(ch_chunks, by: 0)

        GLIMPSE2_SPLITREFERENCE( split_input, ch_map )
        ch_versions = ch_versions.mix( GLIMPSE2_SPLITREFERENCE.out.versions.first() )

        ch_chunks_panel_map = GLIMPSE2_SPLITREFERENCE.out.bin_ref
            .combine( ch_map, by: 0 )
            .map{ metaP, bin_ref, gmap -> [metaP, [], [], bin_ref, [], gmap] }
    } else {
        ch_chunks_panel_map = ch_chunks
            .combine(ch_ref, by:0)
            .combine(ch_map, by:0)
    }

    ch_phase_input = ch_input
        .combine(ch_chunks_panel_map)
        .map{ metaI, bam, bai, bamlist, samples, metaCPM, regionin, regionout, panel, panel_index, gmap ->
            [
                metaI + metaCPM, // combined metadata
                bam, bai, bamlist, samples, // input files
                regionin, regionout, // chunk regions
                panel, panel_index, gmap // panel and map files
            ]
        }

    ch_phase_input.ifEmpty{
        error "BAM_VCF_IMPUTE_GLIMPSE2: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input."
    }

    // Impute with Glimpse2
    GLIMPSE2_PHASE(ch_phase_input, ch_fasta)
    ch_versions = ch_versions.mix( GLIMPSE2_PHASE.out.versions.first() )

    // Index phased file
    BCFTOOLS_INDEX_1(GLIMPSE2_PHASE.out.phased_variants)
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_1.out.versions.first() )

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE2_PHASE.out.phased_variants
        .join( BCFTOOLS_INDEX_1.out.csi )
        .groupTuple()

    GLIMPSE2_LIGATE( ligate_input )
    ch_versions = ch_versions.mix( GLIMPSE2_LIGATE.out.versions.first() )

    BCFTOOLS_INDEX_2( GLIMPSE2_LIGATE.out.merged_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_2.out.versions.first() )

    // Join imputed and index files
    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants
        .join(BCFTOOLS_INDEX_2.out.tbi.mix(BCFTOOLS_INDEX_2.out.csi))

    emit:
    ch_chunks     = ch_chunks          // channel: [ val(meta), txt ]
    ch_vcf_index  = ch_vcf_index       // channel: [ val(meta), vcf, csi ]

    versions      = ch_versions        // channel: [ versions.yml ]
}
