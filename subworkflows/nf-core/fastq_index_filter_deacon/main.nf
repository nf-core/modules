include { DEACON_INDEX  } from '../../../modules/nf-core/deacon/index/main'
include { DEACON_FILTER } from '../../../modules/nf-core/deacon/filter/main'

workflow FASTQ_INDEX_FILTER_DEACON {

    take:
    ch_fasta_reads // [ val(meta), [ fasta ], [ reads ] ]

    main:

    // Check if fastqs are single-end or paired-end and run Deacon accordingly
    ch_reads = ch_fasta_reads
        .map  { meta, fasta, reads ->
            if (meta.single_end) {
                if (reads instanceof List && reads.size() != 1) {
                    error("Error: Check your meta.single_end value. Single-end reads should contain one file only.")
                }
                return [ meta, fasta, reads ]
            } else {
                if (!(reads instanceof List) || reads.size() != 2) {
                    error("Error: Check your meta.single_end value. Paired-end reads should contain two files; a forward and a reverse.")
                }
                return [ meta, fasta, reads ]
            }
        }

    // Extract unique reference fasta files and create fasta-specific metadata
    // This ensures each unique reference is indexed only once
    ch_unique_fastas = ch_reads
        .map { _meta, fasta, _reads -> fasta }
        .unique()
        .map { fasta ->
            def fasta_meta = [ id: fasta.baseName ]
            [ fasta_meta, fasta ]
        }

    // Index unique FASTA files only
    DEACON_INDEX ( ch_unique_fastas )

    // Match indexes back to original samples
    ch_indexes = DEACON_INDEX.out.index
        .map { fasta_meta, index -> [ fasta_meta.id, index ] }
    ch_reads_with_index = ch_reads
        .map { meta, fasta, reads ->
            [ fasta.baseName, meta, reads ]
        }
        .combine(ch_indexes, by: 0)
        .map { _fasta_id, meta, reads, index ->
            [ meta, index, reads ]
        }

    // Filter reads using the matched index
    DEACON_FILTER(ch_reads_with_index)

    // TODO: optionally create output channel with indexes and their original sample-level metadata,
    // this preserves the original behaviour of the workflow
    // ch_index_with_meta = ch_reads_with_index
    //     .map { meta, index, _reads -> [ meta, index ] }

    emit:
    index          = DEACON_INDEX.out.index           // channel: [ val(meta), [ index ] ]  // TODO: optional emit ch_index_with_meta instead
    fastq_filtered = DEACON_FILTER.out.fastq_filtered // channel: [ val(meta), [ fastq ] ]
    summary        = DEACON_FILTER.out.log            // channel: [ val(meta), [ log ] ]
}
