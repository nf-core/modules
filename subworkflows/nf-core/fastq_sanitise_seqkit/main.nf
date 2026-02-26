include { SEQKIT_SANA } from '../../../modules/nf-core/seqkit/sana/main'
include { SEQKIT_PAIR } from '../../../modules/nf-core/seqkit/pair/main'

workflow FASTQ_SANITISE_SEQKIT {

    take:
    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = channel.empty()

    // Add strandness information to meta
    ch_reads_with_strandness = ch_reads
        // seqkit/sana can only receive one file at a time
        .flatMap  { meta, reads ->
            if (meta.single_end) {
                if (reads instanceof List && reads.size() != 1) {
                    error("Error: Check your meta.single_end value. Single-end reads should contain one file only.")
                }
                return [[ meta + [strandness: 'single'], reads ]]
            } else {
                if (!(reads instanceof List) || reads.size() != 2) {
                    error("Error: Check your meta.single_end value. Paired-end reads should contain two files; a forward and a reverse.")
                }
                return [
                    [ meta + [strandness: 'R1'], reads[0] ],
                    [ meta + [strandness: 'R2'], reads[1] ]
                ]
            }
        }

    SEQKIT_SANA( ch_reads_with_strandness )
    ch_versions = ch_versions.mix(SEQKIT_SANA.out.versions.first())

    ch_sanitized_reads = SEQKIT_SANA.out.reads
        .map { meta, fastq ->
            // Remove strandness field from meta to merge back together
            def clean_meta = meta.findAll { key, _value -> key != 'strandness' }
            return [ clean_meta, fastq ]
        }
        .groupTuple(by: 0)
        .branch {
            meta, fastq ->
                single_end: meta.single_end
                    return [ meta, fastq ]
                paired_end: !meta.single_end
                    return [ meta, fastq ]
        }

    SEQKIT_PAIR ( ch_sanitized_reads.paired_end )
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions.first())

    ch_reads = ch_sanitized_reads.single_end.mix(SEQKIT_PAIR.out.reads, SEQKIT_PAIR.out.unpaired_reads)

    emit:
    reads    = ch_reads    // channel: [ val(meta), [ fastq ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
