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
    ANGSD_DOCOUNTS ( ch_bam.combine(channel.of([[]])) )

    ANGSD_CONTAMINATION (
        ANGSD_DOCOUNTS.out.icounts,
        ch_hapmap_file
    )

    emit:
    txt      = ANGSD_CONTAMINATION.out.txt // channel: [ val(meta), [ txt ] ]
}
