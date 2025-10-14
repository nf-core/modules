include { SEQKIT_SANA } from '../../../modules/nf-core/seqkit/sana/main'
include { SEQKIT_PAIR } from '../../../modules/nf-core/seqkit/pair/main'

workflow FASTQ_SEQKIT_SANA_PAIR {

    take:
    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    SEQKIT_SANA( ch_reads.transpose() ) // seqkit/sana can only receive one file at a time
    ch_versions = ch_versions.mix(SEQKIT_SANA.out.versions.first())

    ch_sanitized_reads = SEQKIT_SANA.out.reads
        .groupTuple(by: 0)
        .branch {
            meta, fastq ->
                single_end: meta.single_end
                    return [ meta, fastq ]
                paired_end: !meta.single_end
                    return [ meta, fastq ]
        }

    ch_sanitized_paired_reads = ch_sanitized_reads.paired_end
        .map { meta, reads ->
            def sorted_reads = reads.sort { it.size() } // sorting by file size, to keep snapshot order same
            def renamed = sorted_reads.indexed().collect { i, file_path ->
                def base = file(file_path).getSimpleName()
                def ext  = file(file_path).getName() - file(file_path).getSimpleName()
                def new_name = "${base}_${i + 1}${ext}"
                file(file_path).copyTo(new_name)
            }
            [ meta, renamed ]
        }

    SEQKIT_PAIR ( ch_sanitized_paired_reads )
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions.first())

    ch_reads = ch_sanitized_reads.single_end.mix(SEQKIT_PAIR.out.reads, SEQKIT_PAIR.out.unpaired_reads)

    emit:
    reads    = ch_reads    // channel: [ val(meta), [ fastq ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
