/*
 * Peak calling workflow with HOMER
 */
include { HOMER_MAKETAGDIRECTORY                                } from '../../../modules/nf-core/homer/maketagdirectory/main'
include { HOMER_MAKEUCSCFILE                                    } from '../../../modules/nf-core/homer/makeucscfile/main'
include { HOMER_FINDPEAKS                                       } from '../../../modules/nf-core/homer/findpeaks/main'
include { HOMER_MERGEPEAKS                                      } from '../../../modules/nf-core/homer/mergepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_INDIVIDUAL } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MERGED     } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_POS2BED                                         } from '../../../modules/nf-core/homer/pos2bed/main'
include { UNZIP                                                 } from '../../../modules/nf-core/unzip/main'

workflow HOMER_PEAKCALLING {
    take:
    bam                      // channel: [ val(meta), path(bam) ]
    fasta                    // channel: path(fasta)
    gtf                      // channel: path(gtf) or empty channel. If empty channel, no annotation step
    control                  // channel: [ val(meta), path(control_bam) ] or empty channel
    uniqmap                  // channel: path(uniqmap) or empty channel
    merge_peaks              // val: boolean - whether to merge peaks across samples
    annotate_individual      // val: boolean - whether to annotate individual peak files

    main:
    if(!gtf && annotate_individual){
        log.warn "Individual peak annotation requested but no GTF file provided - skipping annotation"
    }
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

    //
    // Create control tag directory if provided
    //
    ch_control_tagdir = Channel.empty()
    if (control) {
        HOMER_MAKETAGDIRECTORY(
            control,
            fasta
        )
        ch_control_tagdir = HOMER_MAKETAGDIRECTORY.out.tagdir
        ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY.out.versions.first())
    }

    //
    // Create tag directories for all samples
    //
    HOMER_MAKETAGDIRECTORY(
        bam,
        fasta
    )
    ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY.out.versions.first())

    //
    // Creating UCSC Visualization Files
    //
    HOMER_MAKEUCSCFILE(HOMER_MAKETAGDIRECTORY.out.tagdir)
    ch_versions = ch_versions.mix(HOMER_MAKEUCSCFILE.out.versions)

    //
    // Call peaks with or without control
    //
    ch_findpeaks_input = ch_control_tagdir
        ? HOMER_MAKETAGDIRECTORY.out.tagdir.combine(ch_control_tagdir.map { it[1] })
        : HOMER_MAKETAGDIRECTORY.out.tagdir.map { meta, tagdir -> [meta, tagdir, []] }

    HOMER_FINDPEAKS(
        ch_findpeaks_input,
        ch_uniqmap
    )
    ch_versions = ch_versions.mix(HOMER_FINDPEAKS.out.versions.first())

    //
    // Convert peaks to BED format
    //
    HOMER_POS2BED(HOMER_FINDPEAKS.out.txt)
    ch_versions = ch_versions.mix(HOMER_POS2BED.out.versions.first())

    //
    // Annotate individual peaks if requested
    //
    ch_annotated_individual = Channel.empty()
    if (annotate_individual && gtf) {
        HOMER_ANNOTATEPEAKS_INDIVIDUAL(
            HOMER_FINDPEAKS.out.txt,
            fasta,
            gtf
        )
        ch_annotated_individual = HOMER_ANNOTATEPEAKS_INDIVIDUAL.out.txt
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_INDIVIDUAL.out.versions.first())
    }

    //
    // Merge peaks and annotate if requested
    //
    ch_merged_peaks = Channel.empty()
    ch_annotated_merged = Channel.empty()

    if (merge_peaks) {
        ch_all_peaks = HOMER_FINDPEAKS.out.txt
            .map { _meta, txt -> txt }
            .collect()
            .map { peaks -> [[id: 'merged_peaks'], peaks] }

        HOMER_MERGEPEAKS(ch_all_peaks)
        ch_merged_peaks = HOMER_MERGEPEAKS.out.txt
        ch_versions = ch_versions.mix(HOMER_MERGEPEAKS.out.versions)

        if(gtf){
            HOMER_ANNOTATEPEAKS_MERGED(
                HOMER_MERGEPEAKS.out.txt,
                fasta,
                gtf
            )
            ch_annotated_merged = HOMER_ANNOTATEPEAKS_MERGED.out.txt
            ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MERGED.out.versions)
        }
    }

    emit:
    tagdir               = HOMER_MAKETAGDIRECTORY.out.tagdir    // channel: [ val(meta), path(tagdir) ]
    bed_graph            = HOMER_MAKEUCSCFILE.out.bedGraph      // channel: [ val(meta), [ tag_dir/*ucsc.bedGraph.gz ] ]
    peaks                = HOMER_FINDPEAKS.out.txt              // channel: [ val(meta), path(peaks.txt) ]
    bed                  = HOMER_POS2BED.out.bed                // channel: [ val(meta), path(peaks.bed) ]
    merged_peaks         = ch_merged_peaks                      // channel: [ val(meta), path(merged.txt) ] or empty
    annotated_individual = ch_annotated_individual              // channel: [ val(meta), path(annotated.txt) ] or empty
    annotated_merged     = ch_annotated_merged                  // channel: [ val(meta), path(annotated.txt) ] or empty
    versions             = ch_versions                          // channel: [ versions.yml ]
}
