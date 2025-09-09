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

    ch_versions = Channel.empty()

    ch_uniqmap = Channel.empty()

    if (!uniqmap) {
        ch_uniqmap = []
    }
    else if (uniqmap.endsWith('.zip')) {
        ch_uniqmap = UNZIP([[:], uniqmap]).unzipped_archive.map { it[1] }
        ch_versions = ch_versions.mix(UNZIP.out.versions)
    }
    else {
        ch_uniqmap = uniqmap
    }

    /*
    * Create a Tag Directory From The GRO-Seq experiment
    */
    HOMER_MAKETAGDIRECTORY(bam, fasta)
    ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY.out.versions)

    /*
    * Creating UCSC Visualization Files
    */
    HOMER_MAKEUCSCFILE(HOMER_MAKETAGDIRECTORY.out.tagdir)
    ch_versions = ch_versions.mix(HOMER_MAKEUCSCFILE.out.versions)

    /*
    * Find transcripts directly from GRO-Seq
    */
    HOMER_FINDPEAKS(HOMER_MAKETAGDIRECTORY.out.tagdir, ch_uniqmap)
    ch_versions = ch_versions.mix(HOMER_FINDPEAKS.out.versions)

    /*
    * Convert peak file to bed file
    */
    HOMER_POS2BED(HOMER_FINDPEAKS.out.txt)
    ch_versions = ch_versions.mix(HOMER_POS2BED.out.versions)

    emit:
    tagdir    = HOMER_MAKETAGDIRECTORY.out.tagdir // channel: [ val(meta), [ tagdir ] ]
    bed_graph = HOMER_MAKEUCSCFILE.out.bedGraph // channel: [ val(meta), [ tag_dir/*ucsc.bedGraph.gz ] ]
    peaks     = HOMER_FINDPEAKS.out.txt // channel: [ val(meta), [ *peaks.txt ] ]
    bed       = HOMER_POS2BED.out.bed // channel: [ val(meta), [ *peaks.txt ] ]
    versions  = ch_versions // channel: [ versions.yml ]
}
