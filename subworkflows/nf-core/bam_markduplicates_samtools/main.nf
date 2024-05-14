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
    ch_fasta // channel: /path/to/fasta


    main:
    ch_versions = Channel.empty()


    SAMTOOLS_COLLATE ( ch_bam, ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATE.out.versions.first())

    SAMTOOLS_FIXMATE ( SAMTOOLS_COLLATE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions.first())

    SAMTOOLS_SORT ( SAMTOOLS_FIXMATE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_MARKDUP ( SAMTOOLS_SORT.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions.first())


    emit:
    bam      = SAMTOOLS_MARKDUP.out.bam        // channel: [ val(meta), [ bam ] ]
    versions = ch_versions                     // channel: [ versions.yml ]

}
