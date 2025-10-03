//
// GENERATE BED FILE OF GAPS AND LENGTH IN REFERENCE
//

include { SEQTK_CUTN    } from '../../../modules/nf-core/seqtk/cutn/main'
include { GAWK          } from '../../../modules/nf-core/gawk/main'

workflow GAP_FINDER {
    take:
    ch_reference     // Channel [ val(meta), path(fasta) ]

    main:
    ch_versions     = Channel.empty()


    //
    // MODULE: GENERATES A GAP SUMMARY FILE
    //
    SEQTK_CUTN (
        ch_reference
    )
    ch_versions     = ch_versions.mix( SEQTK_CUTN.out.versions )


    //
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE - INPUT FOR PRETEXT MODULE
    //
    GAWK (
        SEQTK_CUTN.out.bed,
        [],
        false
    )
    ch_versions     = ch_versions.mix( GAWK.out.versions )


    emit:
    gap_file        = GAWK.out.output
    versions        = ch_versions
}
