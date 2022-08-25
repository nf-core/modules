include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../../modules/custom/sratoolsncbisettings/main'
include { SRATOOLS_PREFETCH           } from '../../../modules/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP        } from '../../../modules/sratools/fasterqdump/main'

/**
 * Download FASTQ sequencing reads from the NCBI's Sequence Read Archive (SRA).
 */
workflow SRAFASTQ {
    take:
    sra_ids  // channel: [ val(meta), val(id) ]

    main:

    ch_versions = Channel.empty()

    //
    // Detect existing NCBI user settings or create new ones.
    //
    CUSTOM_SRATOOLSNCBISETTINGS()
    def settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings
    ch_versions = ch_versions.mix( CUSTOM_SRATOOLSNCBISETTINGS.out.versions )

    //
    // Prefetch sequencing reads in SRA format.
    //
    SRATOOLS_PREFETCH ( sra_ids, settings )
    ch_versions = ch_versions.mix( SRATOOLS_PREFETCH.out.versions.first() )

    //
    // Convert the SRA format into one or more compressed FASTQ files.
    //
    SRATOOLS_FASTERQDUMP ( SRATOOLS_PREFETCH.out.sra, settings )
    ch_versions = ch_versions.mix( SRATOOLS_FASTERQDUMP.out.versions.first() )

    emit:
    reads    = SRATOOLS_FASTERQDUMP.out.reads  // channel: [ val(meta), [ reads ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
