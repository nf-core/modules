/*
 * Identify transcripts with homer
 */

include { HOMER_MAKETAGDIRECTORY } from '../../../../modules/nf-core/homer/maketagdirectory/main'
include { HOMER_MAKEUCSCFILE     } from '../../../../modules/nf-core/homer/makeucscfile/main'
include { HOMER_FINDPEAKS        } from '../../../../modules/nf-core/homer/findpeaks/main'
include { HOMER_POS2BED          } from '../../../../modules/nf-core/homer/pos2bed/main'

workflow HOMER_GROSEQ {
    take:
    ch_bam   // channel: [ val(meta), path(reads) ]
    ch_fasta // channel: [ path(fasta)]   (file: /path/to/bwa/index/)

    main:

    ch_versions = Channel.empty()

    /*
    * Create a Tag Directory From The GRO-Seq experiment
    */
    HOMER_MAKETAGDIRECTORY ( ch_bam, ch_fasta )
    ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY.out.versions.first())

    /*
    * Creating UCSC Visualization Files
    */
    HOMER_MAKEUCSCFILE ( HOMER_MAKETAGDIRECTORY.out.tagdir )
    ch_versions = ch_versions.mix(HOMER_MAKEUCSCFILE.out.versions.first())

    /*
    * Find transcripts directly from GRO-Seq
    */
    HOMER_FINDPEAKS ( HOMER_MAKETAGDIRECTORY.out.tagdir )
    ch_versions = ch_versions.mix(HOMER_FINDPEAKS.out.versions.first())

    /*
    * Convert peak file to bed file
    */
    HOMER_POS2BED ( HOMER_FINDPEAKS.out.txt )
    ch_versions = ch_versions.mix(HOMER_POS2BED.out.versions.first())

    emit:
    tagdir             = HOMER_MAKETAGDIRECTORY.out.tagdir // channel: [ val(meta), path(tagdir) ]
    bed_graph          = HOMER_MAKEUCSCFILE.out.bedGraph   // channel: [ val(meta), path(bed_graph) ]
    peaks              = HOMER_FINDPEAKS.out.txt           // channel: [ val(meta), path(peaks) ]
    bed                = HOMER_POS2BED.out.bed             // channel: [ val(meta), path(bed) ]

    versions = ch_versions                                 // channel: [ path(versions.yml) ]
}
