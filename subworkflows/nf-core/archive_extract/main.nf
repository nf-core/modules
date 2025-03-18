include { GUNZIP } from '../../../modules/nf-core/gunzip'
include { UNTAR  } from '../../../modules/nf-core/untar'
include { UNZIP  } from '../../../modules/nf-core/unzip'

workflow ARCHIVE_EXTRACT {
    take:
    archive // Channel: [[meta], archive]

    main:
    versions = Channel.empty()

    archive_to_extract = archive.branch { _meta, archive_ ->
        tar: archive_.toString().endsWith('.tar.gz')
        gz: archive_.toString().endsWith('.gz')
        zip: archive_.toString().endsWith('.zip')
        non_assigned: true
    }

    // This is a confidence check
    not_extracted = archive_to_extract.non_assigned
    not_extracted.view { _meta, archive_ -> log.warn("Archive not in the expected format: " + archive_) }

    // Extract any archive with a recognized extension
    GUNZIP(archive_to_extract.gz)
    UNTAR(archive_to_extract.tar)
    UNZIP(archive_to_extract.zip)

    extracted = Channel
        .empty()
        .mix(
            GUNZIP.out.gunzip,
            UNTAR.out.untar,
            UNZIP.out.unzipped_archive,
        )

    versions = versions.mix(GUNZIP.out.versions)
    versions = versions.mix(UNTAR.out.versions)
    versions = versions.mix(UNZIP.out.versions)

    emit:
    extracted     // channel: [ meta, extracted_archive ]
    not_extracted // channel: [ meta, not_extracted_archive ]
    versions      // channel: [ versions.yml ]
}
