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
        .multiMap { meta, vcf, csi, pedigree, region ->
            metadata: [ meta.id, meta]
            region  : [ meta.id, region]
        }

    // Create the File in bed format and use the meta id for the file name
    ch_merged_region = ch_region.region
        .collectFile { metaid, region -> ["${metaid}.bed", region.replace(":","\t").replace("-","\t")] }
        .map { file -> [file.baseName, file] }

    // Link back the meta map with the file
    ch_region_file = ch_region.metadata
        .join(ch_merged_region, failOnMismatch:true, failOnDuplicate:true)
        .map { mid, meta, region_file -> [meta, region_file]}

    BEDTOOLS_MAKEWINDOWS(ch_region_file)
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())

    ch_chunk_output = BEDTOOLS_MAKEWINDOWS.out.bed
        .splitCsv(header: ['Chr', 'Start', 'End'], sep: "\t", skip: 0)
        .map { meta, it -> [meta, it["Chr"]+":"+it["Start"]+"-"+it["End"]]}

    // Count the number of chunks
    ch_chunks_number = BEDTOOLS_MAKEWINDOWS.out.bed
        .map { meta, bed -> [meta, bed.countLines().intValue()]}

    ch_phase_input = ch_vcf
        .map { meta, vcf, index, pedigree, region ->
            [meta, vcf, index, pedigree] }
        .combine(ch_chunk_output, by:0)
        .map { meta, vcf, index, pedigree, chunk ->
                [meta + [id: "${meta.id}_${chunk.replace(":","-")}"], // The meta.id field need to be modified to be unique for each chunk
                vcf, index, pedigree, chunk]}

    SHAPEIT5_PHASECOMMON ( ch_phase_input,
                            ch_ref,
                            ch_scaffold,
                            ch_map )
    ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions.first())

    VCF_INDEX1(SHAPEIT5_PHASECOMMON.out.phased_variant)
    ch_versions = ch_versions.mix(VCF_INDEX1.out.versions.first())

    ch_ligate_input = SHAPEIT5_PHASECOMMON.out.phased_variant
        .join(VCF_INDEX1.out.csi, failOnMismatch:true, failOnDuplicate:true)
        .view()
        .map{ meta, vcf, csi ->
            newmeta = meta + [id: meta.id.split("_")[0..-2].join("_")]
            [newmeta, vcf, csi]}.view()
        .combine(ch_chunks_number, by:0)
        .map{meta, vcf, csi, chunks_num ->
            [groupKey(meta, chunks_num), vcf, csi]}
        .groupTuple()
        .map{ meta, vcf, csi ->
                [ meta,
                vcf.sort { a, b ->
                    def aStart = a.getName().split('-')[-2].toInteger()
                    def bStart = b.getName().split('-')[-2].toInteger()
                    aStart <=> bStart},
                csi]}

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

