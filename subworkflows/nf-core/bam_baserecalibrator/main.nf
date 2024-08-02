//
// PREPARE RECALIBRATION
//

include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'

workflow BAM_BASERECALIBRATOR {
    take:
    ch_cram            // channel: [mandatory] [ meta, cram_markduplicates, crai ]
    ch_dict            // channel: [mandatory] [ dict ]
    ch_fasta           // channel: [mandatory] [ fasta ]
    ch_fasta_fai       // channel: [mandatory] [ fasta_fai ]
    ch_intervals       // channel: [mandatory] [ intervals, num_intervals ] (or [ [], 0 ] if no intervals)
    ch_known_sites     // channel: [optional]  [ known_sites ]
    ch_known_sites_tbi // channel: [optional]  [ known_sites_tbi ]

    main:
    ch_versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    ch_cram_intervals = ch_cram.combine(ch_intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ] }

    // RUN BASERECALIBRATOR
    GATK4_BASERECALIBRATOR(
        ch_cram_intervals,
        ch_fasta.map{ it -> [ it ] },
        ch_fasta_fai.map{ it -> [ it ] },
        ch_dict.map{ it -> [ it ] },
        ch_known_sites,
        ch_known_sites_tbi
    )

    // Figuring out if there is one or more table(s) from the same sample
    ch_table_to_merge = GATK4_BASERECALIBRATOR.out.table.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple().branch{
        // Use meta.num_intervals to asses number of intervals
        single:   it[0].num_intervals <= 1
        multiple: it[0].num_intervals > 1
    }

    // Only when using intervals
    GATK4_GATHERBQSRREPORTS( ch_table_to_merge.multiple )

    // Mix intervals and no_intervals channels together
    ch_table_bqsr = GATK4_GATHERBQSRREPORTS.out.table.mix(ch_table_to_merge.single.map{ meta, table -> [ meta, table[0] ] })
        // Remove no longer necessary field: num_intervals
        .map{ meta, table -> [ meta - meta.subMap('num_intervals'), table ] }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)
    ch_versions = ch_versions.mix(GATK4_GATHERBQSRREPORTS.out.versions)

    emit:
    table_bqsr = ch_table_bqsr // channel: [ meta, table ]

    versions   = ch_versions   // channel: [ versions.yml ]
}
