include { BEDTOOLS_MAKEWINDOWS              } from '../../../modules/nf-core/bedtools/makewindows/main.nf'
include { SHAPEIT5_PHASECOMMON              } from '../../../modules/nf-core/shapeit5/phasecommon/main'
include { SHAPEIT5_LIGATE                   } from '../../../modules/nf-core/shapeit5/ligate/main'
include { BCFTOOLS_INDEX as VCF_INDEX1      } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX2      } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_PHASE_SHAPEIT5 {

    take:
    ch_vcf        // channel (mandatory): [ val(meta), path(vcf), path(csi), path(pedigree), val(region) ]
    ch_ref        // channel (optional) : [ val(meta), path(ref), path(csi) ]
    ch_scaffold   // channel (optional) : [ val(meta), path(scaffold), path(csi) ]
    ch_map        // channel (optional) : [ val(meta), path(map)]

    main:

    ch_versions = Channel.empty()

    // It is needed to generate a file containing the region to phase in a Chr \tab Start \tab End format
    // The meta map needing to be conserved the following steps a required

    // Keep the meta map and the region in two separated channel but keed id field to link them back
    ch_region = ch_vcf
        .multiMap { m, vcf, csi, pedigree, region ->
            metadata: [ m.id, m]
            region  : [ m.id, region]
        }

    // Create the File in bed format and use the meta id for the file name
    ch_merged_region = ch_region.region
        .collectFile { mid, region -> [ "${mid}.bed", region.replace(":","\t").replace("-","\t") ] }
        .map { file -> [ file.baseName, file ] }

    // Link back the meta map with the file
    ch_region_file = ch_region.metadata
        .join(ch_merged_region, failOnMismatch:true, failOnDuplicate:true)
        .map { mid, meta, region_file -> [meta, region_file] }

    BEDTOOLS_MAKEWINDOWS(ch_region_file)
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())

    ch_chunk_output = BEDTOOLS_MAKEWINDOWS.out.bed
        .splitCsv(header: ['Chr', 'Start', 'End'], sep: "\t", skip: 0)
        .map { meta, it -> [meta, it["Chr"]+":"+it["Start"]+"-"+it["End"]] }

    ch_phase_input = ch_vcf
        .map { m, vcf, index, pedigree, region ->
            [m, vcf, index, pedigree] }
        .combine(ch_chunk_output, by:0)

    SHAPEIT5_PHASECOMMON ( ch_phase_input,
                            ch_ref,
                            ch_scaffold,
                            ch_map )
    ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions.first())

    VCF_INDEX1(SHAPEIT5_PHASECOMMON.out.phased_variant)
    ch_versions = ch_versions.mix(VCF_INDEX1.out.versions.first())

    ch_ligate_input = SHAPEIT5_PHASECOMMON.output.phased_variant
        .map{meta, vcf -> [meta, vcf]}
        .groupTuple()
        .combine(VCF_INDEX1.out.csi
            .map{meta, vcf -> [meta, vcf]}
            .groupTuple(),
            by:0)

    SHAPEIT5_LIGATE(ch_ligate_input)
    ch_versions = ch_versions.mix(SHAPEIT5_LIGATE.out.versions.first())

    VCF_INDEX2(SHAPEIT5_LIGATE.out.merged_variants)
    ch_versions = ch_versions.mix(VCF_INDEX2.out.versions.first())

    emit:
    bed                 = BEDTOOLS_MAKEWINDOWS.out.bed           // channel: [ val(meta), bed ]
    variants_phased     = SHAPEIT5_LIGATE.out.merged_variants    // channel: [ val(meta), vcf ]
    variants_index      = VCF_INDEX2.out.csi                     // channel: [ val(meta), csi ]
    versions            = ch_versions                            // channel: [ versions.yml ]
}

