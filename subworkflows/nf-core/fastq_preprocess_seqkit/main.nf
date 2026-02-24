include { FASTQ_SANITISE_SEQKIT } from '../fastq_sanitise_seqkit/main'
include { SEQKIT_SEQ            } from '../../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_REPLACE        } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQKIT_RMDUP          } from '../../../modules/nf-core/seqkit/rmdup/main'

workflow FASTQ_PREPROCESS_SEQKIT {

    take:
    ch_reads               // channel: [ val(meta), [ fastq ] ]
    skip_seqkit_sana_pair  // boolean
    skip_seqkit_seq        // boolean
    skip_seqkit_replace    // boolean
    skip_seqkit_rmdup      // boolean

    main:
    ch_versions = channel.empty()

    if (!skip_seqkit_sana_pair) {
        FASTQ_SANITISE_SEQKIT( ch_reads )
        ch_reads    = FASTQ_SANITISE_SEQKIT.out.reads
        ch_versions = ch_versions.mix(FASTQ_SANITISE_SEQKIT.out.versions)
    }

    // Split paired-end reads and add strandedness to meta
    ch_reads_split = ch_reads
        .flatMap { meta, reads ->
            if (meta.single_end) {
                if (reads instanceof List && reads.size() != 1) {
                    error("Error: Check your meta.single_end value. Single-end reads should contain one file only.")
                }
                return [[ meta + [strandness: 'single'], reads ]]
            } else {
                if (!(reads instanceof List) || reads.size() != 2) {
                    error("Error: Check your meta.single_end value. Paired-end data should have exactly 2 files.")
                }
                return [
                    [ meta + [strandness: 'R1'], reads[0] ],
                    [ meta + [strandness: 'R2'], reads[1] ]
                ]
            }
        }

    if (!skip_seqkit_seq) {
        SEQKIT_SEQ( ch_reads_split )
        ch_reads_split = SEQKIT_SEQ.out.fastx
        ch_versions    = ch_versions.mix(SEQKIT_SEQ.out.versions.first())
    }

    if (!skip_seqkit_replace) {
        SEQKIT_REPLACE( ch_reads_split )
        ch_reads_split = SEQKIT_REPLACE.out.fastx
    }

    if (!skip_seqkit_rmdup) {
        SEQKIT_RMDUP( ch_reads_split )
        ch_reads_split = SEQKIT_RMDUP.out.fastx
        ch_versions    = ch_versions.mix(SEQKIT_RMDUP.out.versions.first())
    }

    ch_reads = ch_reads_split
        .map { meta, fastq ->
            // Remove strandness field from meta to merge back together
            def clean_meta = meta.findAll { key, _value -> key != 'strandness' }
            return [ clean_meta, fastq ]
        }
        .groupTuple(by: 0)
        .map { meta, files ->
            if (meta.single_end) {
                return [ meta, files[0] ]
            } else {
                def sorted_files = files.flatten().sort { index -> index.name }
                return [ meta, sorted_files ]
            }
        }

    emit:
    reads    = ch_reads     // channel: [ val(meta), [ fastq ] ]
    versions = ch_versions  // channel: [ versions.yml ]
}
