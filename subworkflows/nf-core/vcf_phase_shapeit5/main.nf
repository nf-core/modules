include { GLIMPSE2_CHUNK                          } from '../../../modules/nf-core/glimpse2/chunk'
include { SHAPEIT5_PHASECOMMON                    } from '../../../modules/nf-core/shapeit5/phasecommon'
include { SHAPEIT5_LIGATE                         } from '../../../modules/nf-core/shapeit5/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index'

workflow VCF_PHASE_SHAPEIT5 {

    take:
    ch_vcf            // channel (mandatory) : [ [id, chr], vcf, index, pedigree ]
    ch_chunks         // channel (optional)  : [ [chr], regionout ]
    ch_ref            // channel (optional)  : [ [panelid, chr], vcf, index ]
    ch_scaffold       // channel (optional)  : [ [scaffoldid, chr], vcf, index ]
    ch_map            // channel (optional)  : [ [chr], region, map]
    chunk             // val     (mandatory) : boolean to activate/deactivate chunking step
    chunk_model       // val     (mandatory) : model to used for GLIMPSE2_chunk
    array_common_meta // array   (mandatory) : list of meta keys to perform the joins on

    main:

    ch_vcf_map = ch_vcf
        .map { meta, vcf, index, pedigree -> [
            meta.subMap(array_common_meta), meta, vcf, index, pedigree
        ]}
        .combine(
            ch_map
                .map { meta, region, gmap -> [
                    meta.subMap(array_common_meta), meta, region, gmap
                ]},
            by: 0
        )

    if ( chunk == true ){
        // Error if pre-defined chunks are provided when chunking is activated
        ch_chunks
            .filter { _meta, regionout -> regionout.size() > 0 }
            .subscribe {
                error "ERROR: Cannot provide pre-defined chunks (regionin) when chunk=true. Please either set chunk=false to use provided chunks, or remove input chunks to enable automatic chunking."
            }

        GLIMPSE2_CHUNK ( ch_vcf_map.map{
            _metaCommon, metaV, vcf, index, _pedigree, metaM, region, gmap -> [
                metaM + metaV, vcf, index, region, gmap
            ]
        }, chunk_model )

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

    ch_chunks_counts = ch_chunks
        .groupTuple()
        .map { meta, regionouts ->
            [meta, regionouts.size()]
        }

    ch_ref_scaffold_chunks = ch_ref
        .map { meta, vcf, index -> [
            meta.subMap(array_common_meta), meta, vcf, index
        ]}
        .combine(
            ch_scaffold
                .map { meta, vcf, index -> [
                    meta.subMap(array_common_meta), meta, vcf, index
                ]},
            by:0
        )
        .combine(
            ch_chunks
                .combine(ch_chunks_counts, by: 0)
                .map { meta, regionbuf, region_size -> [
                    meta.subMap(array_common_meta), meta, regionbuf, region_size
                ]},
            by:0
        )

    // Make channel with all parameters
    ch_parameters = ch_vcf_map
        .combine(ch_ref_scaffold_chunks, by: 0)

    ch_parameters.ifEmpty{
        error "ERROR: join operation resulted in an empty channel. Please provide a valid ch_map, ch_ref, ch_scaffold and ch_chunks channel as input (same meta map)."
    }

    // Rearrange channel for phasing
    ch_phase_input = ch_parameters
        .map{
            _metaCommon, metaV, vcf, index, pedigree, _metaM, _regionout, gmap, metaR, ref_vcf, ref_index, metaS, scaffold_vcf, scaffold_index, metaC, regionbuf, region_size ->
            def chr = regionbuf.tokenize(':')[0]
            def region = regionbuf.tokenize(':')[1]
            def start = region.tokenize('-')[0]
            def end = region.tokenize('-')[1]
            def paddedStart = String.format('%010d', start as long)
            def paddedEnd = String.format('%010d', end as long)
            def regionoutPadded = "${chr}:${paddedStart}-${paddedEnd}"
            [
                metaR + metaS + metaC + metaV + ["regionout": regionbuf, "regionoutPadded": regionoutPadded, "regionSize": region_size],
                vcf, index,
                pedigree,
                regionbuf,
                ref_vcf, ref_index,
                scaffold_vcf, scaffold_index,
                gmap
            ]
        }

    SHAPEIT5_PHASECOMMON (ch_phase_input)

    BCFTOOLS_INDEX_PHASE(SHAPEIT5_PHASECOMMON.out.phased_variant)

    ch_ligate_input = SHAPEIT5_PHASECOMMON.out.phased_variant
        .join(
            BCFTOOLS_INDEX_PHASE.out.index,
            failOnMismatch:true, failOnDuplicate:true
        )
        .map { meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionout', 'regionoutPadded', 'regionSize']
            [
                groupKey(meta.subMap(keysToKeep), meta.regionSize),
                vcf, index
            ]
        }
        .groupTuple()
        .map { groupKeyObj, vcf, index ->
            // Extract the actual meta from the groupKey
            def meta = groupKeyObj.getGroupTarget()
            [meta, vcf, index]
        }
        .branch { meta, vcf, index ->
            one: vcf.size() == 1
                return [meta, vcf.get(0), index.get(0)]
            more: vcf.size() > 1
                return [meta, vcf, index]
        }

    SHAPEIT5_LIGATE(ch_ligate_input.more, '')

    BCFTOOLS_INDEX_LIGATE(SHAPEIT5_LIGATE.out.merged_variants)

    ch_vcf_index = ch_ligate_input.one
        .mix(SHAPEIT5_LIGATE.out.merged_variants
            .join(
                BCFTOOLS_INDEX_LIGATE.out.index,
                failOnMismatch:true, failOnDuplicate:true
            )
        )

    emit:
    chunks    = ch_chunks    // channel: [ [id, chr], regionout]
    vcf_index = ch_vcf_index // channel: [ [id, chr], vcf, index ]
}
