/*
 * Peak calling workflow with HOMER
 */
include { HOMER_MAKETAGDIRECTORY as HOMER_MAKETAGDIRECTORY_EXP      } from '../../../modules/nf-core/homer/maketagdirectory/main'
include { HOMER_MAKETAGDIRECTORY as HOMER_MAKETAGDIRECTORY_CNTRL    } from'../../../modules/nf-core/homer/maketagdirectory/main'
include { HOMER_MAKEUCSCFILE                                        } from '../../../modules/nf-core/homer/makeucscfile/main'
include { HOMER_FINDPEAKS                                           } from '../../../modules/nf-core/homer/findpeaks/main'
include { HOMER_MERGEPEAKS                                          } from '../../../modules/nf-core/homer/mergepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_INDIVIDUAL     } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MERGED         } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_QUANTIFYPEAKS                } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_POS2BED                                             } from '../../../modules/nf-core/homer/pos2bed/main'
include { UNZIP                                                     } from '../../../modules/nf-core/unzip/main'

workflow HOMER_PEAKCALLING {
    take:
    style                    // val: one of 'factor', 'histone', 'groseq', 'tss', 'dnase', 'super', 'mC'
    bam                      // channel: [ val(meta), path(bam) ] or []
    tagdir                   // channel: [ val(meta), path(tagdir) ] or [] - if provided, skip makeTagDirectory
    fasta                    // channel: path(fasta)
    gtf                      // channel: path(gtf) or []
    control_bam              // channel: [ val(meta), path(control_bam) ] or []
    control_tagdir           // channel: [ val(meta), path(control_tagdir) ] or []
    uniqmap                  // channel: path(uniqmap) or []
    merge_peaks              // val: boolean - whether to merge peaks across samples
    annotate_individual      // val: boolean - whether to annotate individual peak files
    quantify_peaks           // val: boolean - whether to quantify reads in merged peaks (creates count matrix)
    make_bedgraph            // val: boolean - whether to create bedGraph files for UCSC visualization

    main:
    ch_versions = Channel.empty()
    ch_uniqmap = Channel.empty()

    // Handle uniqmap
    if (uniqmap) {
        if (uniqmap.toString().endsWith('.zip')) {
            ch_uniqmap = UNZIP([[:], uniqmap]).unzipped_archive.map { it[1] }
            ch_versions = ch_versions.mix(UNZIP.out.versions)
        }
        else {
            ch_uniqmap = uniqmap
        }
    }
    else {
        ch_uniqmap = []
    }

    //
    // Create or use existing tag directories
    //
    ch_sample_tagdir = Channel.empty()
    ch_control_tagdir = Channel.empty()

    // Create or use sample tag directories
    if (tagdir) {
        ch_sample_tagdir = tagdir
    } else {
        HOMER_MAKETAGDIRECTORY_EXP(bam, fasta)
        ch_sample_tagdir = HOMER_MAKETAGDIRECTORY_EXP.out.tagdir
        ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY_EXP.out.versions.first())
    }

    // Handle control separately
    if (control_tagdir) {
        ch_control_tagdir = control_tagdir
    } else if (control_bam) {
        HOMER_MAKETAGDIRECTORY_CNTRL(control_bam, fasta)
        ch_control_tagdir = HOMER_MAKETAGDIRECTORY_CNTRL.out.tagdir
        ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY_CNTRL.out.versions.first())
    }

    //
    // Creating UCSC Visualization Files
    //
    if ( make_bedgraph ){
        HOMER_MAKEUCSCFILE(ch_sample_tagdir)
        ch_versions = ch_versions.mix(HOMER_MAKEUCSCFILE.out.versions.first())
        ch_bedgraph = HOMER_MAKEUCSCFILE.out.bedGraph

    } else {
        ch_bedgraph = Channel.empty()
    }

    //
    // Call peaks with or without control
    //
    HOMER_FINDPEAKS(
        style,
        ch_sample_tagdir,
        ch_control_tagdir,
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
    // Merge peaks, annotate, and quantify if requested
    //
    ch_merged_txt = Channel.empty()
    ch_annotated_merged = Channel.empty()
    ch_count_matrix = Channel.empty()

    if (merge_peaks) {
        // Collect all peak files for merging
        ch_all_peaks = HOMER_FINDPEAKS.out.txt
            .map { _meta, txt -> txt }
            .collect()
            .map { peaks -> [[id: 'merged'], peaks] }

        HOMER_MERGEPEAKS(ch_all_peaks)
        ch_merged_txt = HOMER_MERGEPEAKS.out.txt
        ch_versions = ch_versions.mix(HOMER_MERGEPEAKS.out.versions)

        // Annotate merged peaks (genomic annotation only)
        if (gtf) {
            HOMER_ANNOTATEPEAKS_MERGED(
                HOMER_MERGEPEAKS.out.txt,
                fasta,
                gtf
            )
            ch_annotated_merged = HOMER_ANNOTATEPEAKS_MERGED.out.txt
            ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MERGED.out.versions.first())
        }

        // Quantify reads in merged peaks across all samples (count matrix)
        if (quantify_peaks && gtf) {
            // Create input for quantification: merged peaks + list of all tag directories
            ch_quantify_input = HOMER_MERGEPEAKS.out.txt
                .combine(
                    ch_sample_tagdir
                        .map { it -> it[1] }
                        .collect()
                        .map { tagdirs -> [[id: 'quantification'], tagdirs] }
                )
                .map { peak_meta, peaks, _tagdir_meta, tagdirs ->
                    [peak_meta, peaks, tagdirs]
                }

            HOMER_QUANTIFYPEAKS(
                ch_quantify_input.map { meta, peaks, _tagdirs -> [meta, peaks] },
                fasta,
                gtf
            )
            ch_count_matrix = HOMER_QUANTIFYPEAKS.out.txt
            ch_versions = ch_versions.mix(HOMER_QUANTIFYPEAKS.out.versions.first())
        }
    }

    emit:
    tagdir               = ch_sample_tagdir                 // channel: [ val(meta), path(tagdir) ]
    bedgraph             = ch_bedgraph                      // channel: [ val(meta), path(bedGraph.gz) ]
    txt                  = HOMER_FINDPEAKS.out.txt          // channel: [ val(meta), path(_{regions,peaks,tss,superEnhancers,...}.txt) ]
    gtf                  = HOMER_FINDPEAKS.out.gtf          // channel: [ val(meta), path(.gtf) ] . this is only for style = 'groseq'
    bed                  = HOMER_POS2BED.out.bed            // channel: [ val(meta), path(_{regions,peaks,tss,superEnhancers,...}.bed) ]
    merged_txt           = ch_merged_txt                    // channel: [ val(meta), path(merged.txt) ] or empty
    annotated_individual = ch_annotated_individual          // channel: [ val(meta), path(annotated.txt) ] or empty
    annotated_merged     = ch_annotated_merged              // channel: [ val(meta), path(annotated.txt) ] or empty
    count_matrix         = ch_count_matrix                  // channel: [ val(meta), path(count_matrix.txt) ] or empty
    versions             = ch_versions                      // channel: [ versions.yml ]
}
