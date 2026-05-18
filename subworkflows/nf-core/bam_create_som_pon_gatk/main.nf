//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//

include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/nf-core/gatk4/createsomaticpanelofnormals'
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_MERGEMUTECTSTATS            } from '../../../modules/nf-core/gatk4/mergemutectstats'
include { GATK4_MERGEVCFS                   } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MUTECT2                     } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_CREATE_SOM_PON_GATK {
    take:
    ch_mutect2_in // channel: [ val(meta), path(input), path(input_index) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai), path(gzi) ]
    ch_dict // channel: [ val(meta), path(dict) ]
    val_pon_norm // string:  name for panel of normals
    ch_gendb_intervals // channel: [ path(interval_file) ]
    ch_intervals // channel: [ path(intervals), val(num_intervals) ] or [ [], 0 ] if no intervals

    main:
    //
    // Combine input and intervals for scatter and gather strategy
    // Move num_intervals to meta map for GATK4_MUTECT2
    //
    ch_mutect2_intervals = ch_mutect2_in
        .combine(ch_intervals)
        .map { meta, input, input_index, intervals_, num_intervals ->
            [meta + [num_intervals: num_intervals], input, input_index, intervals_]
        }

    //
    // Perform variant calling for each sample using mutect2 module in panel of normals mode
    //
    GATK4_MUTECT2(
        ch_mutect2_intervals,
        ch_fasta,
        ch_fai,
        ch_dict,
        [],
        [],
        [],
        [],
        [],
        [],
    )

    // Branch outputs based on whether intervals were used
    vcf_branch = GATK4_MUTECT2.out.vcf.branch { meta, _vcf ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    tbi_branch = GATK4_MUTECT2.out.tbi.branch { meta, _tbi ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    stats_branch = GATK4_MUTECT2.out.stats.branch { meta, _stats ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Only when using intervals: group outputs by sample and merge
    vcf_to_merge = vcf_branch.intervals.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple()
    stats_to_merge = stats_branch.intervals.map { meta, stats -> [groupKey(meta, meta.num_intervals), stats] }.groupTuple()

    GATK4_MERGEVCFS(vcf_to_merge, ch_dict)
    GATK4_MERGEMUTECTSTATS(stats_to_merge)

    // Mix merged and non-interval outputs, remove num_intervals from meta
    ch_vcf = channel.empty()
        .mix(GATK4_MERGEVCFS.out.vcf, vcf_branch.no_intervals)
        .map { meta, vcf -> [meta - meta.subMap('num_intervals'), vcf] }
        .collect { _meta, vcf -> [vcf] }
        .toList()

    ch_tbi = channel.empty()
        .mix(GATK4_MERGEVCFS.out.tbi, tbi_branch.no_intervals)
        .map { meta, tbi -> [meta - meta.subMap('num_intervals'), tbi] }
        .collect { _meta, tbi -> [tbi] }
        .toList()

    //
    // Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport
    //
    ch_dict_gendb = ch_dict.map { _meta, dict -> [dict] }.toList()

    ch_gendb_input = channel.of([id: val_pon_norm])
        .combine(ch_vcf)
        .combine(ch_tbi)
        .combine(ch_gendb_intervals)
        .combine(ch_dict_gendb)
        .map { meta, vcf, tbi, interval, dict -> [meta, vcf, tbi, interval, [], dict] }

    GATK4_GENOMICSDBIMPORT(ch_gendb_input, false, false, false)

    //
    // Panel of normals made from genomicsdb workspace using createsomaticpanelofnormals
    //
    GATK4_CREATESOMATICPANELOFNORMALS(GATK4_GENOMICSDBIMPORT.out.genomicsdb, ch_fasta, ch_fai.map { meta, fai, _gzi -> [meta, fai] }, ch_dict)

    emit:
    genomicsdb = GATK4_GENOMICSDBIMPORT.out.genomicsdb // channel: [ val(meta), path(genomicsdb) ]
    pon_vcf    = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf // channel: [ val(meta), path(vcf) ]
    pon_index  = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi // channel: [ val(meta), path(tbi) ]
}
