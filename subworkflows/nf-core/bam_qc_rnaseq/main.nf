//
// Run post-alignment QC tools on RNA-seq BAM files
//

include { DUPRADAR                        } from '../../../modules/nf-core/dupradar/main'
include { PRESEQ_LCEXTRAP                 } from '../../../modules/nf-core/preseq/lcextrap/main'
include { QUALIMAP_RNASEQ                 } from '../../../modules/nf-core/qualimap/rnaseq/main'
include { SUBREAD_FEATURECOUNTS           } from '../../../modules/nf-core/subread/featurecounts/main'
include { CUSTOM_MULTIQCCUSTOMBIOTYPE     } from '../../../modules/nf-core/custom/multiqccustombiotype/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QUALIMAP } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_RSEQC                       } from '../bam_rseqc/main'

workflow BAM_QC_RNASEQ {

    take:
    ch_bam_bai              // channel: [ val(meta), path(bam), path(bai) ]
    ch_gtf                  // channel: [ val(meta), path(gtf) ]
    ch_gene_bed             // channel: path(bed)
    ch_fasta_fai            // channel: [ val(meta), path(fasta), path(fai) ]
    ch_biotypes_header      // channel: [ val(meta), path(biotypes_header.txt) ]
    skip_preseq             // val(bool)
    skip_biotype_qc         // val(bool)
    skip_qualimap           // val(bool)
    skip_dupradar           // val(bool)
    skip_rseqc              // val(bool)
    biotype                 // val(string) - e.g. "gene_type" or "gene_biotype"
    rseqc_modules           // val(list)   - e.g. ['bam_stat', 'infer_experiment', ...]

    main:
    ch_multiqc_files       = channel.empty()
    ch_inferexperiment_txt = channel.empty()

    ch_genome_bam = ch_bam_bai.map { meta, bam, _bai -> [ meta, bam ] }

    // Initialise individual output channels (populated conditionally below)
    ch_preseq_lc_extrap       = channel.empty()
    ch_preseq_log             = channel.empty()
    ch_biotype_tsv            = channel.empty()
    ch_biotype_rrna           = channel.empty()
    ch_featurecounts_counts   = channel.empty()
    ch_featurecounts_summary  = channel.empty()
    ch_qualimap_results       = channel.empty()
    ch_dupradar_scatter2d     = channel.empty()
    ch_dupradar_boxplot       = channel.empty()
    ch_dupradar_hist          = channel.empty()
    ch_dupradar_dupmatrix     = channel.empty()
    ch_dupradar_intercept_slope = channel.empty()
    ch_dupradar_multiqc       = channel.empty()
    ch_bamstat_txt              = channel.empty()
    ch_innerdistance_all        = channel.empty()
    ch_innerdistance_distance   = channel.empty()
    ch_innerdistance_freq       = channel.empty()
    ch_innerdistance_mean       = channel.empty()
    ch_innerdistance_pdf        = channel.empty()
    ch_innerdistance_rscript    = channel.empty()
    ch_junctionannotation_all          = channel.empty()
    ch_junctionannotation_bed          = channel.empty()
    ch_junctionannotation_interact_bed = channel.empty()
    ch_junctionannotation_xls          = channel.empty()
    ch_junctionannotation_pdf          = channel.empty()
    ch_junctionannotation_events_pdf   = channel.empty()
    ch_junctionannotation_rscript      = channel.empty()
    ch_junctionannotation_log          = channel.empty()
    ch_junctionsaturation_all     = channel.empty()
    ch_junctionsaturation_pdf     = channel.empty()
    ch_junctionsaturation_rscript = channel.empty()
    ch_readdistribution_txt     = channel.empty()
    ch_readduplication_all      = channel.empty()
    ch_readduplication_seq_xls  = channel.empty()
    ch_readduplication_pos_xls  = channel.empty()
    ch_readduplication_pdf      = channel.empty()
    ch_readduplication_rscript  = channel.empty()
    ch_tin_txt                  = channel.empty()

    //
    // MODULE: Run Preseq
    //
    if (!skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_preseq_lc_extrap = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_preseq_log       = PRESEQ_LCEXTRAP.out.log
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap)
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    if (!skip_biotype_qc && biotype) {

        SUBREAD_FEATURECOUNTS (
            ch_genome_bam
                .combine(ch_gtf.map { _meta, gtf -> gtf })
        )

        ch_featurecounts_counts  = SUBREAD_FEATURECOUNTS.out.counts
        ch_featurecounts_summary = SUBREAD_FEATURECOUNTS.out.summary

        CUSTOM_MULTIQCCUSTOMBIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header
        )
        ch_biotype_tsv  = CUSTOM_MULTIQCCUSTOMBIOTYPE.out.tsv
        ch_biotype_rrna = CUSTOM_MULTIQCCUSTOMBIOTYPE.out.rrna
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_MULTIQCCUSTOMBIOTYPE.out.tsv)
    }

    //
    // MODULE: Qualimap
    //
    if (!skip_qualimap) {
        // Sort BAM by name for qualimap
        // Requires ext.args = '-n' to be set by the caller
        SAMTOOLS_SORT_QUALIMAP (
            ch_genome_bam,
            ch_fasta_fai,
            ''
        )

        QUALIMAP_RNASEQ (
            SAMTOOLS_SORT_QUALIMAP.out.bam,
            ch_gtf
        )
        ch_qualimap_results = QUALIMAP_RNASEQ.out.results
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results)
    }

    //
    // MODULE: dupRadar
    //
    if (!skip_dupradar) {
        DUPRADAR (
            ch_genome_bam,
            ch_gtf
        )
        ch_dupradar_scatter2d     = DUPRADAR.out.scatter2d
        ch_dupradar_boxplot       = DUPRADAR.out.boxplot
        ch_dupradar_hist          = DUPRADAR.out.hist
        ch_dupradar_dupmatrix     = DUPRADAR.out.dupmatrix
        ch_dupradar_intercept_slope = DUPRADAR.out.intercept_slope
        ch_dupradar_multiqc       = DUPRADAR.out.multiqc
        ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc)
    }

    //
    // SUBWORKFLOW: RSeQC
    //
    if (!skip_rseqc && rseqc_modules.size() > 0) {
        BAM_RSEQC (
            ch_bam_bai.map { meta, bam, bai -> [ meta, [ bam, bai ] ] },
            ch_gene_bed,
            rseqc_modules
        )
        ch_bamstat_txt              = BAM_RSEQC.out.bamstat_txt
        ch_innerdistance_all        = BAM_RSEQC.out.innerdistance_all
        ch_innerdistance_distance   = BAM_RSEQC.out.innerdistance_distance
        ch_innerdistance_freq       = BAM_RSEQC.out.innerdistance_freq
        ch_innerdistance_mean       = BAM_RSEQC.out.innerdistance_mean
        ch_innerdistance_pdf        = BAM_RSEQC.out.innerdistance_pdf
        ch_innerdistance_rscript    = BAM_RSEQC.out.innerdistance_rscript
        ch_inferexperiment_txt      = BAM_RSEQC.out.inferexperiment_txt
        ch_junctionannotation_all          = BAM_RSEQC.out.junctionannotation_all
        ch_junctionannotation_bed          = BAM_RSEQC.out.junctionannotation_bed
        ch_junctionannotation_interact_bed = BAM_RSEQC.out.junctionannotation_interact_bed
        ch_junctionannotation_xls          = BAM_RSEQC.out.junctionannotation_xls
        ch_junctionannotation_pdf          = BAM_RSEQC.out.junctionannotation_pdf
        ch_junctionannotation_events_pdf   = BAM_RSEQC.out.junctionannotation_events_pdf
        ch_junctionannotation_rscript      = BAM_RSEQC.out.junctionannotation_rscript
        ch_junctionannotation_log          = BAM_RSEQC.out.junctionannotation_log
        ch_junctionsaturation_all     = BAM_RSEQC.out.junctionsaturation_all
        ch_junctionsaturation_pdf     = BAM_RSEQC.out.junctionsaturation_pdf
        ch_junctionsaturation_rscript = BAM_RSEQC.out.junctionsaturation_rscript
        ch_readdistribution_txt     = BAM_RSEQC.out.readdistribution_txt
        ch_readduplication_all      = BAM_RSEQC.out.readduplication_all
        ch_readduplication_seq_xls  = BAM_RSEQC.out.readduplication_seq_xls
        ch_readduplication_pos_xls  = BAM_RSEQC.out.readduplication_pos_xls
        ch_readduplication_pdf      = BAM_RSEQC.out.readduplication_pdf
        ch_readduplication_rscript  = BAM_RSEQC.out.readduplication_rscript
        ch_tin_txt                  = BAM_RSEQC.out.tin_txt

        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.bamstat_txt)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.inferexperiment_txt)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.innerdistance_freq)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionannotation_log)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionsaturation_rscript)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readdistribution_txt)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readduplication_pos_xls)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.tin_txt)
    }

    emit:
    // Aggregated
    multiqc_files                       = ch_multiqc_files                               // channel: [ val(meta), path(files) ]

    // Preseq
    preseq_lc_extrap                    = ch_preseq_lc_extrap                            // channel: [ val(meta), path(txt) ]
    preseq_log                          = ch_preseq_log                                  // channel: [ val(meta), path(log) ]

    // Biotype QC
    featurecounts_counts                = ch_featurecounts_counts                        // channel: [ val(meta), path(txt) ]
    featurecounts_summary               = ch_featurecounts_summary                       // channel: [ val(meta), path(txt) ]
    biotype_tsv                         = ch_biotype_tsv                                 // channel: [ val(meta), path(tsv) ]
    biotype_rrna                        = ch_biotype_rrna                                // channel: [ val(meta), path(tsv) ]

    // Qualimap
    qualimap_results                    = ch_qualimap_results                            // channel: [ val(meta), path(dir) ]

    // dupRadar
    dupradar_scatter2d                  = ch_dupradar_scatter2d                          // channel: [ val(meta), path(pdf) ]
    dupradar_boxplot                    = ch_dupradar_boxplot                            // channel: [ val(meta), path(pdf) ]
    dupradar_hist                       = ch_dupradar_hist                               // channel: [ val(meta), path(pdf) ]
    dupradar_dupmatrix                  = ch_dupradar_dupmatrix                          // channel: [ val(meta), path(txt) ]
    dupradar_intercept_slope            = ch_dupradar_intercept_slope                    // channel: [ val(meta), path(txt) ]
    dupradar_multiqc                    = ch_dupradar_multiqc                            // channel: [ val(meta), path(txt) ]

    // RSeQC (populated via ch_ variables set inside the conditional block)
    inferexperiment_txt                 = ch_inferexperiment_txt                          // channel: [ val(meta), path(txt) ]
    bamstat_txt                         = ch_bamstat_txt                                  // channel: [ val(meta), path(txt) ]
    innerdistance_all                   = ch_innerdistance_all                            // channel: [ val(meta), path(txt/pdf/r) ]
    innerdistance_distance              = ch_innerdistance_distance                       // channel: [ val(meta), path(txt) ]
    innerdistance_freq                  = ch_innerdistance_freq                           // channel: [ val(meta), path(txt) ]
    innerdistance_mean                  = ch_innerdistance_mean                           // channel: [ val(meta), path(txt) ]
    innerdistance_pdf                   = ch_innerdistance_pdf                            // channel: [ val(meta), path(pdf) ]
    innerdistance_rscript               = ch_innerdistance_rscript                        // channel: [ val(meta), path(r) ]
    junctionannotation_all              = ch_junctionannotation_all                       // channel: [ val(meta), path(bed/xls/pdf/r/log) ]
    junctionannotation_bed              = ch_junctionannotation_bed                       // channel: [ val(meta), path(bed) ]
    junctionannotation_interact_bed     = ch_junctionannotation_interact_bed              // channel: [ val(meta), path(bed) ]
    junctionannotation_xls              = ch_junctionannotation_xls                       // channel: [ val(meta), path(xls) ]
    junctionannotation_pdf              = ch_junctionannotation_pdf                       // channel: [ val(meta), path(pdf) ]
    junctionannotation_events_pdf       = ch_junctionannotation_events_pdf                // channel: [ val(meta), path(pdf) ]
    junctionannotation_rscript          = ch_junctionannotation_rscript                   // channel: [ val(meta), path(r) ]
    junctionannotation_log              = ch_junctionannotation_log                       // channel: [ val(meta), path(log) ]
    junctionsaturation_all              = ch_junctionsaturation_all                       // channel: [ val(meta), path(pdf/r) ]
    junctionsaturation_pdf              = ch_junctionsaturation_pdf                       // channel: [ val(meta), path(pdf) ]
    junctionsaturation_rscript          = ch_junctionsaturation_rscript                   // channel: [ val(meta), path(r) ]
    readdistribution_txt                = ch_readdistribution_txt                         // channel: [ val(meta), path(txt) ]
    readduplication_all                 = ch_readduplication_all                          // channel: [ val(meta), path(xls/pdf/r) ]
    readduplication_seq_xls             = ch_readduplication_seq_xls                      // channel: [ val(meta), path(xls) ]
    readduplication_pos_xls             = ch_readduplication_pos_xls                      // channel: [ val(meta), path(xls) ]
    readduplication_pdf                 = ch_readduplication_pdf                          // channel: [ val(meta), path(pdf) ]
    readduplication_rscript             = ch_readduplication_rscript                       // channel: [ val(meta), path(r) ]
    tin_txt                             = ch_tin_txt                                      // channel: [ val(meta), path(txt) ]

}
