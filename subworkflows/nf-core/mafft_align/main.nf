include { MAFFT_ALIGN as MAFFT_ALIGN_MODULE } from '../../../modules/nf-core/mafft/align/main'

workflow MAFFT_ALIGN {

    take:
    ch_fasta // channel: [ val(meta), fasta ]

    main:

    MAFFT_ALIGN_MODULE ( ch_fasta, [[], []], [[], []], [[], []], [[], []], [[], []], [] )

    emit:
    alignment = MAFFT_ALIGN_MODULE.out.fas      // channel: [ val(meta), *.fas ]
}
