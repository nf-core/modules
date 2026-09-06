/*
 * Identify transcripts with homer
 */

include { UNZIP                  } from '../../../modules/nf-core/unzip/main'

include { HOMER_MAKETAGDIRECTORY } from '../../../modules/nf-core/homer/maketagdirectory/main'
include { HOMER_MAKEUCSCFILE     } from '../../../modules/nf-core/homer/makeucscfile/main'
include { HOMER_FINDPEAKS        } from '../../../modules/nf-core/homer/findpeaks/main'
include { HOMER_POS2BED          } from '../../../modules/nf-core/homer/pos2bed/main'

workflow HOMER_GROSEQ {
    take:
    bam     // channel: [ val(meta), [ reads ] ]
    fasta   //    file: /path/to/bwa/index/
    uniqmap

    main:

    ch_uniqmap = channel.empty()

    if (!uniqmap) {
        ch_uniqmap = []
    }
    else if (uniqmap.endsWith('.zip')) {
        ch_uniqmap = UNZIP([[:], uniqmap]).unzipped_archive.map { index -> index[1] }
    }
    else {
        ch_uniqmap = uniqmap
    }

    /*
    * Create a Tag Directory From The GRO-Seq experiment
    */
    HOMER_MAKETAGDIRECTORY(bam, fasta)

    /*
    * Creating UCSC Visualization Files
    */
    HOMER_MAKEUCSCFILE(HOMER_MAKETAGDIRECTORY.out.tagdir)

    /*
    * Find transcripts directly from GRO-Seq
    */
    HOMER_FINDPEAKS(HOMER_MAKETAGDIRECTORY.out.tagdir, ch_uniqmap)

    /*
    * Convert peak file to bed file
    */
    HOMER_POS2BED(HOMER_FINDPEAKS.out.txt)

    emit:
    tagdir    = HOMER_MAKETAGDIRECTORY.out.tagdir // channel: [ val(meta), [ tagdir ] ]
    bed_graph = HOMER_MAKEUCSCFILE.out.bedGraph // channel: [ val(meta), [ tag_dir/*ucsc.bedGraph.gz ] ]
    peaks     = HOMER_FINDPEAKS.out.txt // channel: [ val(meta), [ *peaks.txt ] ]
    bed       = HOMER_POS2BED.out.bed // channel: [ val(meta), [ *peaks.txt ] ]
}
