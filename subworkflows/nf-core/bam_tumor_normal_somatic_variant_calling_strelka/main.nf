include { GATK4_MERGEVCFS as MERGE_STRELKA_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA_SNVS   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { STRELKA_SOMATIC                         } from '../../../modules/nf-core/strelka/somatic/main'

workflow BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_STRELKA {
    take:
    ch_cram          // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi ] manta* are optional
    ch_dict          // channel: [optional]  [ meta, dict ]
    ch_fasta         // channel: [mandatory] [ fasta ]
    ch_fasta_fai     // channel: [mandatory] [ fasta_fai ]
    ch_intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals

    main:
    ch_versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    ch_cram_intervals = ch_cram.combine(ch_intervals)
        // Move num_intervals to meta map
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, intervals, intervals_index, num_intervals -> [ meta + [ num_intervals:num_intervals ], normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, intervals, intervals_index ] }

    STRELKA_SOMATIC (
        ch_cram_intervals,
        ch_fasta,
        ch_fasta_fai
    )

    // Figuring out if there is one or more vcf(s) from the same sample
    ch_vcf_indels = STRELKA_SOMATIC.out.vcf_indels.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    ch_vcf_snvs = STRELKA_SOMATIC.out.vcf_snvs.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    ch_vcf_indels_to_merge = ch_vcf_indels.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    ch_vcf_snvs_to_merge = ch_vcf_snvs.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_STRELKA_INDELS( ch_vcf_indels_to_merge, ch_dict )
    MERGE_STRELKA_SNVS( ch_vcf_snvs_to_merge, ch_dict )

    // Mix intervals and no_intervals channels together
    ch_vcf = Channel.empty().mix(MERGE_STRELKA_INDELS.out.vcf, MERGE_STRELKA_SNVS.out.vcf, ch_vcf_indels.no_intervals, ch_vcf_snvs.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'strelka' ], vcf ] }

    ch_versions = ch_versions.mix(MERGE_STRELKA_SNVS.out.versions)
    ch_versions = ch_versions.mix(MERGE_STRELKA_INDELS.out.versions)
    ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions)

    emit:
    vcf      = ch_vcf

    versions = ch_versions
}
