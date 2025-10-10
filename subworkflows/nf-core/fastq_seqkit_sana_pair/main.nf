include { SEQKIT_SANA } from '../../../modules/nf-core/seqkit/sana/main'
include { SEQKIT_PAIR } from '../../../modules/nf-core/seqkit/pair/main'

workflow FASTQ_SEQKIT_SANA_PAIR {

    take:
    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:

    ch_versions = Channel.empty()

    ch_reads = ch_reads.transpose() // seqkit/sana can only receives one file at a time
    
    SEQKIT_SANA( ch_reads )
    ch_versions = ch_versions.mix(SEQKIT_SANA.out.versions.first())

    ch_sanitized_reads = SEQKIT_SANA.out.reads
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
        }

    ch_sanitized_reads_paired = ch_sanitized_reads.paired_end
        .groupTuple(by: 0)
        .map { meta, reads ->
            def renamed = reads.indexed().collect { i, file_path ->
                def base = file(file_path).getSimpleName()
                def ext  = file(file_path).getName() - file(file_path).getSimpleName()
                def new_name = "${base}_${i + 1}${ext}"
                file(file_path).copyTo(new_name)
            }
            [ meta, renamed ]
        }

    SEQKIT_PAIR ( ch_sanitized_reads_paired )
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions.first())

    ch_sanitized_reads.single_end.mix(SEQKIT_PAIR.out.reads, SEQKIT_PAIR.out.unpaired_reads).view()

    emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
