//
// RECALIBRATE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_MERGE as MERGE_BAM  } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_MERGE as MERGE_CRAM } from '../../../modules/nf-core/samtools/merge'
include { GATK4_APPLYBQSR  } from '../../../modules/nf-core/gatk4/applybqsr'

workflow BAM_APPLYBQSR {
    take:
    reads       // channel: [mandatory] [ meta, reads, index, recal ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    dict      // channel: [mandatory] [ dict ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    reads_interval = reads
        .combine(intervals)
        .map { meta, reads_, index, recal, intervals_, num_intervals ->
            [meta + [num_intervals: num_intervals], reads_, index, recal, intervals_]
        }

    output_suffix = reads_interval
        .map { _meta, reads_, _index, _recal, _intervals -> reads_.getExtension() == 'bam' ? 'bam' : 'cram' }

    // RUN APPLYBQSR
    GATK4_APPLYBQSR(
        reads_interval,
        fasta.map { _meta, fasta_file -> [fasta_file] },
        fasta_fai.map { _meta, fai_file -> [fai_file] },
        dict.map { _meta, dict_file -> [dict_file] },
        output_suffix,
    )
    // BAM path — populated when ext.suffix='bam', empty otherwise

    bam_applybqsr = GATK4_APPLYBQSR.out.bam
            .join(GATK4_APPLYBQSR.out.bai, failOnDuplicate: true)
            .branch {
                single: it[0].num_intervals == 1
                multiple: it[0].num_intervals > 1
            }

        // For multiple intervals, gather and merge the recalibrated cram files
        bam_to_merge = bam_applybqsr.multiple
            .map { meta, bam_, _bai -> [groupKey(meta, meta.num_intervals), bam_] }
            .groupTuple()

        MERGE_BAM(
                bam_to_merge,
                [ [ id:'null' ], [], [] ]
            )

        bam_recal = MERGE_BAM.out.bam
                .join(MERGE_BAM.out.index, failOnDuplicate: true, failOnMismatch: true)
                .mix(bam_applybqsr.single)
                .map { meta, bam_, bai -> [meta - meta.subMap('num_intervals'), bam_, bai] }

     // ---- CRAM path (populated when a cram was produced) ----
        cram_applybqsr = GATK4_APPLYBQSR.out.cram
            .join(GATK4_APPLYBQSR.out.bai, failOnDuplicate: true)
            .branch {
                single:   it[0].num_intervals == 1
                multiple: it[0].num_intervals > 1
            }

        cram_to_merge = cram_applybqsr.multiple
            .map { meta, cram_, crai -> [groupKey(meta, meta.num_intervals), cram_, crai] }
            .groupTuple()

        MERGE_CRAM(
            cram_to_merge,
            [ [ id:'null' ], [], [] ]
        )

        cram_recal = MERGE_CRAM.out.cram
            .join(MERGE_CRAM.out.index, failOnDuplicate: true, failOnMismatch: true)
            .mix(cram_applybqsr.single)
            .map { meta, cram_, crai -> [meta - meta.subMap('num_intervals'), cram_, crai] }

        // Unified output — bam or cram depending on what was produced
        recal = bam_recal.mix(cram_recal)

        emit:
        recal    // channel: [ meta, file, index ]
    }
}
