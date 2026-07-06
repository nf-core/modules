include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_MERGE  } from '../../../modules/nf-core/samtools/merge'

workflow BAM_APPLYBQSR {
    take:
    input // channel: [mandatory] [ meta, reads, index, recal ]
    references // channel: [mandatory] [ meta, fasta, fai, dict ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    output_suffix // string: [mandatory] 'bam' or 'cram'

    main:
    reads_interval = input
        .combine(intervals)
        .map { meta, reads_, index, recal, intervals_, num_intervals ->
            [meta + [num_intervals: num_intervals], reads_, index, recal, intervals_]
        }

    // RUN APPLYBQSR
    GATK4_APPLYBQSR(reads_interval, references, output_suffix)

    reads_applybqsr = GATK4_APPLYBQSR.out.bam
        .mix(GATK4_APPLYBQSR.out.cram)
        .combine(GATK4_APPLYBQSR.out.bai, by: 0)
        .branch { meta, _reads, _index ->
            single: meta.num_intervals <= 1
            multiple: meta.num_intervals > 1
        }

    // For multiple intervals, gather and merge the recalibrated cram files
    reads_to_merge = reads_applybqsr.multiple
        .map { meta, reads, index -> [groupKey(meta, meta.num_intervals), reads, index] }
        .groupTuple()

    SAMTOOLS_MERGE(reads_to_merge, references.map { meta, fasta, fai, _dict -> [meta, fasta, fai, []] })

    // Unified output — bam or cram depending on what was produced
    recal_out = SAMTOOLS_MERGE.out.bam
        .mix(SAMTOOLS_MERGE.out.cram)
        .combine(SAMTOOLS_MERGE.out.index, by: 0)
        .mix(reads_applybqsr.single)
        .map { meta, reads, index -> [meta - meta.subMap('num_intervals'), reads, index] }

    emit:
    recal_out = recal_out // channel: [ meta, file, index ]
}
