//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//

include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/nf-core/gatk4/createsomaticpanelofnormals/main'
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_MERGEVCFS                   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MUTECT2                     } from '../../../modules/nf-core/gatk4/mutect2/main'

workflow BAM_CREATE_SOM_PON_GATK {
    take:
    ch_input            // channel: [ val(meta), path(input), path(index) ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fai              // channel: [ val(meta), path(fai) ]
    ch_dict             // channel: [ val(meta), path(dict) ]
    val_pon_norm        // string:  name for panel of normals
    ch_intervals        // channel: [mandatory] [ intervals, num_intervals ]

    main:
    ch_versions = Channel.empty()

    //
    // Perform variant calling for each sample using mutect2 module in panel of normals mode.
    //

    // Combine input and intervals for spread and gather strategy
    ch_input_intervals = ch_input.combine(ch_intervals)
        // Move num_intervals to meta map and reorganize channel for MUTECT2_PAIRED module
        // .map{ meta, input, index, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input, index, intervals ] }
        .map{ meta, input, index, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input, index, [] ] }

    GATK4_MUTECT2(ch_input_intervals, ch_fasta, ch_fai, ch_dict, [], [], [], [])

    // GATK4_MERGEVCFS(GATK4_MUTECT2.out.vcf.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple(), ch_dict)

    // vcf_joint = GATK4_MERGEVCFS.out.vcf.map{ meta, vcf -> [[id: 'joint'], vcf]}.groupTuple()
    // tbi_joint = GATK4_MERGEVCFS.out.tbi.map{ meta, tbi -> [[id: 'joint'], tbi]}.groupTuple()

    vcf_joint = GATK4_MUTECT2.out.vcf.map{ meta, vcf -> [[id: 'joint'], vcf]}.groupTuple()
    tbi_joint = GATK4_MUTECT2.out.tbi.map{ meta, tbi -> [[id: 'joint'], tbi]}.groupTuple()

    ch_genomicsdb_input = vcf_joint.join(tbi_joint).combine(ch_intervals).map{ meta, vcf, tbi, intervals, num_intervals -> [meta, vcf, tbi, intervals, [], [] ]}

    GATK4_GENOMICSDBIMPORT(ch_genomicsdb_input, [], [], [])
    GATK4_CREATESOMATICPANELOFNORMALS(GATK4_GENOMICSDBIMPORT.out.genomicsdb, ch_fasta, ch_fai, ch_dict)

    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())
    // ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS.out.versions.first())

    emit:
    vcf          = GATK4_MUTECT2.out.vcf                     // channel: [ val(meta), path(vcf) ]
    index        = GATK4_MUTECT2.out.tbi                     // channel: [ val(meta), path(tbi) ]
    stats        = GATK4_MUTECT2.out.stats                   // channel: [ val(meta), path(stats) ]
    genomicsdb   = GATK4_GENOMICSDBIMPORT.out.genomicsdb     // channel: [ val(meta), path(genomicsdb) ]
    pon_vcf      = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf // channel: [ val(meta), path(vcf) ]
    pon_index    = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi // channel: [ val(meta), path(tbi) ]

    versions     = ch_versions                               // channel: [ path(versions.yml) ]
}
