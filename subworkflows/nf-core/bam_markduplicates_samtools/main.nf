//
// Collate, fixmate, sort and markdup using Samtools
//

include { SAMTOOLS_COLLATE } from '../../../modules/nf-core/samtools/collate/main'
include { SAMTOOLS_FIXMATE } from '../../../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_SORT    } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MARKDUP } from '../../../modules/nf-core/samtools/markdup/main'


workflow BAM_MARKDUPLICATES_SAMTOOLS {

    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), [ fasta ], [fai] ]

    main:

    SAMTOOLS_COLLATE ( ch_bam, ch_fasta )

    SAMTOOLS_FIXMATE ( SAMTOOLS_COLLATE.out.bam )

    SAMTOOLS_SORT (
        SAMTOOLS_FIXMATE.out.bam,
        ch_fasta.map{meta, fasta, _fai -> [meta, fasta]},
        ''
    )

    SAMTOOLS_MARKDUP ( SAMTOOLS_SORT.out.bam, ch_fasta )

    emit:
    bam      = SAMTOOLS_MARKDUP.out.bam        // channel: [ val(meta), [ bam ] ]

}
