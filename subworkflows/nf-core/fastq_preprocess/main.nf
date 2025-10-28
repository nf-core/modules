include { FASTQ_SANITISE_SEQKIT             } from '../fastq_sanitise_seqkit/main'
include { SEQKIT_SEQ                        } from '../../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_REPLACE                    } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQKIT_RMDUP                      } from '../../../modules/nf-core/seqkit/rmdup/main'

workflow FASTQ_PREPROCESS {

    take:
    ch_reads          // channel: [ val(meta), [ fastq ] ]
    skip_seqkit_sana_pair // boolean
    skip_seqkit_seq       // boolean
    skip_seqkit_replace   // boolean
    skip_seqkit_rmdup     // boolean

    main:
    ch_versions = Channel.empty()

    if (!skip_seqkit_sana_pair) {
        FASTQ_SANITISE_SEQKIT( ch_reads )
        ch_reads        = FASTQ_SANITISE_SEQKIT.out.reads
        ch_versions     = ch_versions.mix(FASTQ_SANITISE_SEQKIT.out.versions.first())
    }

    // Split paired-end reads and add unique IDs
    ch_reads_split = ch_reads
        .flatMap { meta, reads ->
            if (meta.single_end) {
                return [[ meta, reads ]]
            } else {
                if (!(reads instanceof List) || reads.size() != 2) {
                    error("Paired-end data should have exactly 2 files")
                }
                return [
                    [ meta + [read_pair: 'R1', id: "${meta.id}_R1"], reads[0] ],
                    [ meta + [read_pair: 'R2', id: "${meta.id}_R2"], reads[1] ]
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
        ch_versions    = ch_versions.mix(SEQKIT_REPLACE.out.versions.first())
    }

    if (!skip_seqkit_rmdup) {
        SEQKIT_RMDUP( ch_reads_split )
        ch_reads_split = SEQKIT_RMDUP.out.fastx
        ch_versions    = ch_versions.mix(SEQKIT_RMDUP.out.versions.first())
    }

    // Merge back and restore original IDs
    ch_reads = ch_reads_split
        .map { meta, fastq ->
            // Restore original ID and remove read_pair field
            def clean_meta = meta.findAll { key, value -> key != 'read_pair' }
            if (!clean_meta.single_end) {
                clean_meta.id = clean_meta.id.replaceAll(/_R[12]$/, '')
            }
            return [ clean_meta, fastq ]
        }

    ch_reads
        .groupTuple(by: 0)
        .map { meta, files ->
            if (meta.single_end) {
                return [ meta, files[0] ]
            } else {
                def sorted_files = files.flatten().sort()
                return [ meta, sorted_files ]
            }
        }
        .set { ch_reads }

    emit:
    reads    = ch_reads    // channel: [ val(meta), [ fastq ] ]
    versions = ch_versions // channel: [ versions.yml ]

}
