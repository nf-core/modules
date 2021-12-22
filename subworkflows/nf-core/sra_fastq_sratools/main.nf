//
// Download FASTQ sequencing reads from the NCBI's Sequence Read Archive (SRA).
//

include { SRATOOLS_PREFETCH    } from '../../../modules/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP } from '../../../modules/sratools/fasterqdump/main'

workflow SRA_FASTQ_SRATOOLS {
    take:
    sra_ids  // channel: [ val(meta), val(id) ]

    main:

    ch_versions = Channel.empty()

    //
    // Prefetch sequencing reads in SRA format.
    //
    SRATOOLS_PREFETCH ( sra_ids )
    ch_versions = ch_versions.mix(SRATOOLS_PREFETCH.out.versions.first())

    //
    // Convert the SRA format into one or more compressed FASTQ files.
    //
    SRATOOLS_FASTERQDUMP ( SRATOOLS_PREFETCH.out.sra )
    ch_versions = ch_versions.mix(SRATOOLS_FASTERQDUMP.out.versions.first())

    emit:
    reads    = SRATOOLS_FASTERQDUMP.out.reads  // channel: [ val(meta), [ reads ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
