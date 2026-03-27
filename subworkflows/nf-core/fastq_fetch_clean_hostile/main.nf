include { HOSTILE_FETCH } from '../../../modules/nf-core/hostile/fetch/main'
include { HOSTILE_CLEAN } from '../../../modules/nf-core/hostile/clean/main'

workflow FASTQ_FETCH_CLEAN_HOSTILE {

    take:
    ch_reads        // channel: [ val(meta), [ fastq ] ]
    ch_reference    // channel: [ val(reference_name), path(reference_dir) ] (optional)
    index_name      // val (optional)

    main:
    ch_versions = channel.empty()

    if (!index_name && !ch_reference) {
        error "Provide either the reference index name for HOSTILE_FETCH or an existing reference path for HOSTILE_CLEAN."
    }

    if (index_name && ch_reference) {
        error "Cannot provide both the reference index name for fetching and a reference path, please provide only one."
    }

    if (index_name) {
        HOSTILE_FETCH( index_name )
        out_reference = HOSTILE_FETCH.out.reference
        ch_versions = ch_versions.mix(HOSTILE_FETCH.out.versions.first())
    }
    else {
        out_reference = ch_reference
    }

    HOSTILE_CLEAN( ch_reads, out_reference )
    ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions.first())

    emit:
        reference   = out_reference                  // channel: [ val(reference_name), path(reference_dir) ]
        fastq       = HOSTILE_CLEAN.out.fastq       // channel: [ val(meta), [ *.fastq.gz ] ]
        json        = HOSTILE_CLEAN.out.json        // channel: [ val(meta), [ *.json ] ]
        versions    = ch_versions                   // channel: [ versions.yml ]
}
