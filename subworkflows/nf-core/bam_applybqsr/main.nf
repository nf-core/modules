include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_MERGE  } from '../../../modules/nf-core/samtools/merge'

workflow BAM_APPLYBQSR {
    take:
    input // channel: [mandatory] [ meta, reads, index, recal ]
    references // channel: [mandatory] [ meta, fasta, fai, dict ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    output_format // tuple: [mandatory] [ file_type, index_format ] — ['bam','bai'], ['bam','csi'], or ['cram','crai']

    main:
    def file_type = output_format[0]
    def index_fmt = output_format[1]

    def valid_formats = [
        'bam'  : ['bai', 'csi'],
        'cram' : ['crai']
    ]

    if (!valid_formats.containsKey(file_type)) {
        error("BAM_APPLYBQSR: invalid file type '${file_type}'. Must be 'bam' or 'cram'.")
    }
    if (!valid_formats[file_type].contains(index_fmt)) {
        error("BAM_APPLYBQSR: invalid index format '${index_fmt}' for '${file_type}' output. Must be one of: ${valid_formats[file_type]}.")
    }

    reads_interval = input
        .combine(intervals)
        .map { meta, reads_, index, recal, intervals_, num_intervals ->
            [meta + [num_intervals: num_intervals, interval: intervals_.name], reads_, index, recal, intervals_]
        }

    // RUN APPLYBQSR
    GATK4_APPLYBQSR(reads_interval, references, file_type)

    reads_applybqsr = GATK4_APPLYBQSR.out.bam
        .mix(GATK4_APPLYBQSR.out.cram)
        .join(GATK4_APPLYBQSR.out.bai.mix(GATK4_APPLYBQSR.out.crai), remainder: true)
        .branch { meta, _reads, _index ->
            single: meta.num_intervals <= 1
            multiple: meta.num_intervals > 1
        }

    // For multiple intervals, gather and merge the recalibrated cram files
    reads_to_merge = reads_applybqsr.multiple
        .map { meta, reads, index ->
            def merge_meta = meta - meta.subMap('interval')
            [groupKey(merge_meta, merge_meta.num_intervals), reads, index]
        }
        .groupTuple()

    SAMTOOLS_MERGE(reads_to_merge, references.map { meta, fasta, fai, _dict -> [meta, fasta, fai, []] }, index_fmt)

    // Unified output — bam or cram depending on what was produced
    recal_out = SAMTOOLS_MERGE.out.bam
        .mix(SAMTOOLS_MERGE.out.cram)
        .join(SAMTOOLS_MERGE.out.index, failOnMismatch: true)
        .mix(reads_applybqsr.single)
        .map { meta, reads, index -> [meta - meta.subMap('num_intervals', 'interval'), reads, index] }

    emit:
    recal_out = recal_out // channel: [ meta, file, index ]
}
