include { GLIMPSE_CHUNK      } from '../../../modules/nf-core/glimpse/chunk/main'
include { GLIMPSE_PHASE      } from '../../../modules/nf-core/glimpse/phase/main'
include { GLIMPSE_LIGATE     } from '../../../modules/nf-core/glimpse/ligate/main'
include { BCFTOOLS_INDEX     } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_IMPUTE_GLIMPSE {

    take:
    ch_vcf      // channel (mandatory): [ val(meta), path(vcf), path(csi), path(region), path(sample) ]
    ch_ref      // channel (mandatory): [ val(meta), path(vcf), path(csi) ]
    ch_map      // channel  (optional): [ path(map) ]           path to map
    ch_infos    // channel  (optional): [ path(info) ]           sample infos

    main:

    ch_versions = Channel.empty()

    GLIMPSE_CHUNK ( ch_vcf )
    ch_versions = ch_versions.mix(GLIMPSE_CHUNK.out.versions)

    chunk_output = GLIMPSE_CHUNK.out.chunk_chr
                                .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
                                .map { metamap, it -> [metamap, it["RegionIn"], it["RegionOut"]]}
    phase_input = ch_vcf.map{[it[0], it[1], it[2]]}
                            .join(chunk_output)
                            .join(ch_ref)
                            .join(ch_map)
                            .join(ch_infos)

    GLIMPSE_PHASE ( phase_input ) // [meta, vcf, index, regionin, regionout, regionindex, ref, ref_index, map, sample_infos]
    ch_versions = ch_versions.mix(GLIMPSE_PHASE.out.versions.first())

    ligate_input  = GLIMPSE_PHASE.out.phased_variant.groupTuple()

    BCFTOOLS_INDEX ( ligate_input )
    GLIMPSE_LIGATE ( ligate_input.join(BCFTOOLS_INDEX.out.csi.groupTuple()) )

    ch_versions = ch_versions.mix(GLIMPSE_LIGATE.out.versions.first())

    emit:
    chunk_chr        = GLIMPSE_CHUNK.out.chunk_chr           // channel: [ val(meta), path(txt) ]
    merged_variants  = GLIMPSE_LIGATE.out.merged_variants    // channel: [ val(meta), path(bcf) ]
    phased_variants  = GLIMPSE_PHASE.out.phased_variant      // channel: [ val(meta), path(bcf) ]

    versions         = ch_versions                           // channel: [ path(versions.yml) ]
}
