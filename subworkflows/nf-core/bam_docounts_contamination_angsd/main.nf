//
// ANGSD doCounts and contamination estimation on the X-chromosome
//

include { ANGSD_DOCOUNTS      } from '../../../modules/nf-core/angsd/docounts/main'
include { ANGSD_CONTAMINATION } from '../../../modules/nf-core/angsd/contamination/main'

workflow BAM_DOCOUNTS_CONTAMINATION_ANGSD {

    take:
    ch_bam         // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_hapmap_file // channel: [ val(meta), [ hapmap_file ] ]

    main:
    ch_versions = channel.empty()

    ANGSD_DOCOUNTS ( ch_bam.combine(channel.of([[]])) )
    ch_versions = ch_versions.mix(ANGSD_DOCOUNTS.out.versions)

    ANGSD_CONTAMINATION (
        ANGSD_DOCOUNTS.out.icounts,
        ch_hapmap_file
    )
    ch_versions = ch_versions.mix(ANGSD_CONTAMINATION.out.versions)

    emit:
    txt      = ANGSD_CONTAMINATION.out.txt // channel: [ val(meta), [ txt ] ]
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}
