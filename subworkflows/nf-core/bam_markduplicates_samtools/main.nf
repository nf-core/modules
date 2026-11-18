// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

//
// Collate, fixmate, sort and markdup using Samtools
//

include { SAMTOOLS_MULTICOMMAND } from '../../../modules/nf-core/samtools/multicommand/main'


workflow BAM_MARKDUPLICATES_SAMTOOLS {
    take:
    ch_bam // channel: [ val(meta), [ bam ],  [ index ] ]
    ch_fasta_fai // channel: [ val(meta), [ fasta ], [fai] ]

    main:

    SAMTOOLS_MULTICOMMAND(
        ch_bam,
        ch_fasta_fai,
        ["collate", "fixmate", "sort", "markdup"],
    )

    emit:
    bam = SAMTOOLS_MULTICOMMAND.out.bam // channel: [ val(meta), [ bam ] ]
}
