//
// ANGSD doCounts and contamination estimation on the X-chromosome
//

include { ANGSD_DOCOUNTS      } from '../../../modules/nf-core/angsd/docounts/main'
include { ANGSD_CONTAMINATION } from '../../../modules/nf-core/angsd/contamination/main'

workflow BAM_DOCOUNTS_CONTAMINATION_ANGSD {

    take:
    ch_bam         // channel: [ val(meta), [ bam ] ]
    ch_bai         // channel: [ val(meta), [ bai ] ]
    ch_hapmap_file // channel: [ val(meta), [ hapmap_file ] ]

    main:
    ch_versions = Channel.empty()

    ch_input = ch_bam
        .join( ch_bai )
        .map{
            meta, bam, bai ->
            [ meta, bam, bai, [] ]
        }

    ANGSD_DOCOUNTS ( ch_input )
    ch_versions = ch_versions.mix(ANGSD_DOCOUNTS.out.versions.first())
    ch_icounts  = ANGSD_DOCOUNTS.out.icounts

    ANGSD_CONTAMINATION ( ch_icounts, ch_hapmap_file )
    ch_versions = ch_versions.mix(ANGSD_CONTAMINATION.out.versions.first())

    emit:
    txt      = ANGSD_CONTAMINATION.out.txt // channel: [ val(meta), [ txt ] ]
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}
