//
// Download FASTQ sequencing reads from the NCBI's Sequence Read Archive (SRA).
//

params.prefetch_options    = [:]
params.fasterqdump_options = [:]

include { SRATOOLS_PREFETCH    } from '../../../modules/sratools/prefetch/main'    addParams( options: params.prefetch_options    )
include { SRATOOLS_FASTERQDUMP } from '../../../modules/sratools/fasterqdump/main' addParams( options: params.fasterqdump_options )

workflow SRA_FASTQ {
    take:
    sra_ids  // channel: [ val(meta), val(id) ]

    main:

    ch_versions = Channel.empty()

    //
    // Prefetch sequencing reads in SRA format.
    //
    SRATOOLS_PREFETCH ( sra_ids )
    ch_versions = ch_versions.mix( SRATOOLS_PREFETCH.out.versions.first() )

    //
    // Convert the SRA format into one or more compressed FASTQ files.
    //
    SRATOOLS_FASTERQDUMP ( SRATOOLS_PREFETCH.out.sra )
    ch_versions = ch_versions.mix( SRATOOLS_FASTERQDUMP.out.versions.first() )

    emit:
    reads    = SRATOOLS_FASTERQDUMP.out.reads  // channel: [ val(meta), [ reads ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
