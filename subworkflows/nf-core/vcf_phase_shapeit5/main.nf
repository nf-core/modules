include { GLIMPSE2_CHUNK                          } from '../../../modules/nf-core/glimpse2/chunk'
include { SHAPEIT5_PHASECOMMON                    } from '../../../modules/nf-core/shapeit5/phasecommon'
include { SHAPEIT5_LIGATE                         } from '../../../modules/nf-core/shapeit5/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index'

workflow VCF_PHASE_SHAPEIT5 {

    take:
    ch_vcf      // channel (mandatory) : [ [id, chr], vcf, index, pedigree, region ]
    ch_chunks   // channel (optional)  : [ [id, chr], regionout ]
    ch_ref      // channel (optional)  : [ [id, chr], vcf, index ]
    ch_scaffold // channel (optional)  : [ [id, chr], vcf, index ]
    ch_map      // channel (optional)  : [ [id, chr], map]
    chunk       // val     (mandatory) : boolean to activate/deactivate chunking step
    chunk_model // channel (mandatory) : [ model ]

    main:

    if ( chunk == true ){
        // Error if pre-defined chunks are provided when chunking is activated
        ch_chunks
            .filter { _meta, regionout -> regionout.size() > 0 }
            .subscribe {
                error "ERROR: Cannot provide pre-defined chunks (regionin) when chunk=true. Please either set chunk=false to use provided chunks, or remove input chunks to enable automatic chunking."
            }

        // Chunk reference panel
        ch_vcf_map = ch_vcf
            .combine(ch_map, by: 0)
            .map{
                meta, vcf, index, _pedigree, region, gmap -> [
                    meta, vcf, index, region, gmap
                ]
            }

        GLIMPSE2_CHUNK ( ch_vcf_map, chunk_model )

        ch_chunks = GLIMPSE2_CHUNK.out.chunk_chr
            .splitCsv(header: [
                'ID', 'Chr', 'RegionBuf', 'RegionCnk', 'WindowCm',
                'WindowMb', 'NbTotVariants', 'NbComVariants'
            ], sep: "\t", skip: 0)
            .map { meta, rows -> [meta, rows["RegionBuf"]]}
    }

    ch_chunks
        .filter { _meta, regionout -> regionout.size() == 0 }
        .subscribe {
            error "ERROR: ch_chunks channel is empty. Please provide a valid channel or set chunk parameter to true."
        }

    // Make channel with all parameters
    ch_parameters = ch_vcf
        .combine(ch_map, by: 0)
        .combine(ch_ref, by: 0)
        .combine(ch_scaffold, by: 0)
        .combine(ch_chunks, by: 0)

    ch_parameters.ifEmpty{
        error "ERROR: join operation resulted in an empty channel. Please provide a valid ch_map, ch_ref, ch_scaffold and ch_chunks channel as input (same meta map)."
    }

    // Rearrange channel for phasing
    ch_phase_input = ch_parameters
        .map{
            meta, vcf, index, pedigree, _region, gmap, ref_vcf, ref_index, scaffold_vcf, scaffold_index, regionbuf -> [
                meta + ["regionout": regionbuf], vcf, index, pedigree, regionbuf,
                ref_vcf, ref_index, scaffold_vcf, scaffold_index, gmap
            ]
        }

    SHAPEIT5_PHASECOMMON (ch_phase_input)

    BCFTOOLS_INDEX_PHASE(SHAPEIT5_PHASECOMMON.out.phased_variant)

    ch_ligate_input = SHAPEIT5_PHASECOMMON.out.phased_variant
        .join(
            BCFTOOLS_INDEX_PHASE.out.tbi.mix(BCFTOOLS_INDEX_PHASE.out.csi),
            failOnMismatch:true, failOnDuplicate:true
        )
        .map{ meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionout']
            [ meta.subMap(keysToKeep), vcf, index ]
        }
        .groupTuple()

    SHAPEIT5_LIGATE(ch_ligate_input)

    BCFTOOLS_INDEX_LIGATE(SHAPEIT5_LIGATE.out.merged_variants)

    ch_vcf_index = SHAPEIT5_LIGATE.out.merged_variants
        .join(
            BCFTOOLS_INDEX_LIGATE.out.tbi.mix(BCFTOOLS_INDEX_LIGATE.out.csi),
            failOnMismatch:true, failOnDuplicate:true
        )

    emit:
    chunks    = ch_chunks    // channel: [ [id, chr], regionout]
    vcf_index = ch_vcf_index // channel: [ [id, chr], vcf, csi ]
}
