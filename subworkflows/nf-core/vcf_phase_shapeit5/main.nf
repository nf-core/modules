include { BEDTOOLS_MAKEWINDOWS              } from '../../../modules/nf-core/bedtools/makewindows/main.nf'
include { SHAPEIT5_PHASECOMMON              } from '../../../modules/nf-core/shapeit5/phasecommon/main'  
include { SHAPEIT5_LIGATE                   } from '../../../modules/nf-core/shapeit5/ligate/main'
include { BCFTOOLS_INDEX as VCF_INDEX1      } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX2      } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_PHASE_SHAPEIT5 {

    take:
    ch_vcf        // channel (mandatory): [ [id, ref], vcf, csi, pedigree, region ]
    ch_ref        // channel (optional) : [ meta, ref, csi ]
    ch_scaffold   // channel (optional) : [ meta, scaffold, csi ]
    ch_map        // channel (optional) : [ meta, map]

    main:

    ch_versions = Channel.empty()

    // It is needed to generate a file containing the region to phase in a Chr \tab Start \tab End format
    // The meta map needing to be conserved the following steps a required

    // Keep the meta map and the region in two separated channel but keed id field to link them back
    ch_vcf
        .multiMap { m, vcf, csi, pedigree, region ->
            metadata: [m.id, m]
            region: [m.id,region ]
        }
        .set { ch_region }

    // Create the File in bed format and use the meta id for the file name
    ch_region.region
        .collectFile { mid, region -> [ "${mid}.bed", region.replace(":","\t").replace("-","\t") ] }
        .map { file -> [ file.simpleName, file ] }
        .set { ch_merged_region }

    // Link back the meta map with the file
    ch_region.metadata
        .join(ch_merged_region)
        .map {mid, meta, region_file -> [meta, region_file]}
        .set { ch_region_file }
    
    BEDTOOLS_MAKEWINDOWS(ch_region_file)
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())

    chunk_output = BEDTOOLS_MAKEWINDOWS.out.bed
                                .splitCsv(header: ['Chr', 'Start', 'End'], sep: "\t", skip: 0)
                                .map { meta, it -> [meta, it["Chr"]+":"+it["Start"]+"-"+it["End"]]}

    phase_input = ch_vcf.map{m, vcf, index, pedigree, region -> [m, vcf, index, pedigree]}
                    .combine(chunk_output, by:0)

    SHAPEIT5_PHASECOMMON ( phase_input,
                            ch_ref,
                            ch_scaffold,
                            ch_map )
    ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions.first())

    VCF_INDEX1(SHAPEIT5_PHASECOMMON.out.phased_variant)
    ch_versions = ch_versions.mix(VCF_INDEX1.out.versions.first())

    ligate_input = SHAPEIT5_PHASECOMMON.output.phased_variant
                                        .map{meta, vcf -> [meta, vcf]}
                                        .groupTuple()
                                        .combine(VCF_INDEX1.out.csi
                                            .map{meta, vcf -> [meta, vcf]}
                                            .groupTuple(),
                                            by:0)

    SHAPEIT5_LIGATE(ligate_input) 
    ch_versions = ch_versions.mix(SHAPEIT5_LIGATE.out.versions.first())

    VCF_INDEX2(SHAPEIT5_LIGATE.out.merged_variants)
    ch_versions = ch_versions.mix(VCF_INDEX2.out.versions.first())

    emit:
    bed                 = BEDTOOLS_MAKEWINDOWS.out.bed           // channel: [ val(meta), bed ]
    variants_phased     = SHAPEIT5_LIGATE.out.merged_variants    // channel: [ val(meta), vcf ]
    variants_index      = VCF_INDEX2.out.csi                     // channel: [ val(meta), csi ]
    versions            = ch_versions                            // channel: [ versions.yml ]
}

