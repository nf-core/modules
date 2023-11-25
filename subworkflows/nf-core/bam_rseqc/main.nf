//
// Run RSeQC modules
//

include { RSEQC_BAMSTAT            } from '../../../modules/nf-core/rseqc/bamstat/main'
include { RSEQC_INNERDISTANCE      } from '../../../modules/nf-core/rseqc/innerdistance/main'
include { RSEQC_INFEREXPERIMENT    } from '../../../modules/nf-core/rseqc/inferexperiment/main'
include { RSEQC_JUNCTIONANNOTATION } from '../../../modules/nf-core/rseqc/junctionannotation/main'
include { RSEQC_JUNCTIONSATURATION } from '../../../modules/nf-core/rseqc/junctionsaturation/main'
include { RSEQC_READDISTRIBUTION   } from '../../../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_READDUPLICATION    } from '../../../modules/nf-core/rseqc/readduplication/main'
include { RSEQC_TIN                } from '../../../modules/nf-core/rseqc/tin/main'

workflow BAM_RSEQC {
    take:
    ch_bam_bai    // channel: [ val(meta), [ bam, bai ] ]
    ch_bed        //    file: /path/to/genome.bed
    rseqc_modules //    list: rseqc modules to run

    main:

    ch_versions = Channel.empty()

    ch_bam_bai
        .map { [ it[0], it[1] ] }
        .set { ch_bam }

    //
    // Run RSeQC bam_stat.py
    //
    bamstat_txt = Channel.empty()
    if ('bam_stat' in rseqc_modules) {
        RSEQC_BAMSTAT ( ch_bam )
        ch_bamstat = RSEQC_BAMSTAT.out.txt
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions.first())
    }

    //
    // Run RSeQC inner_distance.py
    //
    innerdistance_distance = Channel.empty()
    innerdistance_freq     = Channel.empty()
    innerdistance_mean     = Channel.empty()
    innerdistance_pdf      = Channel.empty()
    innerdistance_rscript  = Channel.empty()
    if ('inner_distance' in rseqc_modules) {
        RSEQC_INNERDISTANCE ( ch_bam, ch_bed )
        innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        innerdistance_freq     = RSEQC_INNERDISTANCE.out.freq
        innerdistance_mean     = RSEQC_INNERDISTANCE.out.mean
        innerdistance_pdf      = RSEQC_INNERDISTANCE.out.pdf
        innerdistance_rscript  = RSEQC_INNERDISTANCE.out.rscript
        ch_inner_distance      = innerdistance_distance.mix(innerdistance_freq, innerdistance_mean, innerdistance_pdf, innerdistance_rscript)
        ch_versions = ch_versions.mix(RSEQC_INNERDISTANCE.out.versions.first())
    }

    //
    // Run RSeQC infer_experiment.py
    //
    inferexperiment_txt = Channel.empty()
    if ('infer_experiment' in rseqc_modules) {
        RSEQC_INFEREXPERIMENT ( ch_bam, ch_bed )
        ch_inferexperiment = RSEQC_INFEREXPERIMENT.out.txt
        ch_versions = ch_versions.mix(RSEQC_INFEREXPERIMENT.out.versions.first())
    }

    //
    // Run RSeQC junction_annotation.py
    //
    junctionannotation_bed          = Channel.empty()
    junctionannotation_interact_bed = Channel.empty()
    junctionannotation_xls          = Channel.empty()
    junctionannotation_pdf          = Channel.empty()
    junctionannotation_events_pdf   = Channel.empty()
    junctionannotation_rscript      = Channel.empty()
    junctionannotation_log          = Channel.empty()
    if ('junction_annotation' in rseqc_modules) {
        RSEQC_JUNCTIONANNOTATION ( ch_bam, ch_bed )
        junctionannotation_bed          = RSEQC_JUNCTIONANNOTATION.out.bed
        junctionannotation_interact_bed = RSEQC_JUNCTIONANNOTATION.out.interact_bed
        junctionannotation_xls          = RSEQC_JUNCTIONANNOTATION.out.xls
        junctionannotation_pdf          = RSEQC_JUNCTIONANNOTATION.out.pdf
        junctionannotation_events_pdf   = RSEQC_JUNCTIONANNOTATION.out.events_pdf
        junctionannotation_rscript      = RSEQC_JUNCTIONANNOTATION.out.rscript
        junctionannotation_log          = RSEQC_JUNCTIONANNOTATION.out.log
        ch_junction_annotation          = junctionannotation_bed.mix(junctionannotation_interact_bed, junctionannotation_xls, junctionannotation_pdf, junctionannotation_events_pdf, junctionannotation_rscript, junctionannotation_log)
        ch_versions = ch_versions.mix(RSEQC_JUNCTIONANNOTATION.out.versions.first())
    }

    //
    // Run RSeQC junction_saturation.py
    //
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()
    if ('junction_saturation' in rseqc_modules) {
        RSEQC_JUNCTIONSATURATION ( ch_bam, ch_bed )
        junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
        junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
        ch_junction_saturation     = junctionsaturation_pdf.mix(junctionsaturation_rscript)
        ch_versions = ch_versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())
    }

    //
    // Run RSeQC read_distribution.py
    //
    readdistribution_txt = Channel.empty()
    if ('read_distribution' in rseqc_modules) {
        RSEQC_READDISTRIBUTION ( ch_bam, ch_bed )
        ch_readdistribution = RSEQC_READDISTRIBUTION.out.txt
        ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())
    }

    //
    // Run RSeQC read_duplication.py
    //
    readduplication_seq_xls = Channel.empty()
    readduplication_pos_xls = Channel.empty()
    readduplication_pdf     = Channel.empty()
    readduplication_rscript = Channel.empty()
    if ('read_duplication' in rseqc_modules) {
        RSEQC_READDUPLICATION ( ch_bam )
        readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
        readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
        readduplication_pdf     = RSEQC_READDUPLICATION.out.pdf
        readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
        ch_read_duplication     = readduplication_seq_xls.mix(readduplication_pos_xls, readduplication_pdf, readduplication_rscript)
        ch_versions = ch_versions.mix(RSEQC_READDUPLICATION.out.versions.first())
    }

    //
    // Run RSeQC tin.py
    //
    tin_txt = Channel.empty()
    if ('tin' in rseqc_modules) {
        RSEQC_TIN ( ch_bam_bai, ch_bed )
        ch_tin     = RSEQC_TIN.out.txt
        ch_versions = ch_versions.mix(RSEQC_TIN.out.versions.first())
    }

    emit:
    ch_bamstat                     // channel: [ val(meta), txt ]

    ch_inner_distance
    innerdistance_distance          // channel: [ val(meta), txt ]
    innerdistance_freq              // channel: [ val(meta), txt ]
    innerdistance_mean              // channel: [ val(meta), txt ]
    innerdistance_pdf               // channel: [ val(meta), pdf ]
    innerdistance_rscript           // channel: [ val(meta), r   ]

    ch_inferexperiment             // channel: [ val(meta), txt ]

    ch_junction_annotation
    junctionannotation_bed          // channel: [ val(meta), bed ]
    junctionannotation_interact_bed // channel: [ val(meta), bed ]
    junctionannotation_xls          // channel: [ val(meta), xls ]
    junctionannotation_pdf          // channel: [ val(meta), pdf ]
    junctionannotation_events_pdf   // channel: [ val(meta), pdf ]
    junctionannotation_rscript      // channel: [ val(meta), r   ]
    junctionannotation_log          // channel: [ val(meta), log ]

    ch_junction_saturation
    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    ch_readdistribution            // channel: [ val(meta), txt ]

    ch_read_duplication
    readduplication_seq_xls         // channel: [ val(meta), xls ]
    readduplication_pos_xls         // channel: [ val(meta), xls ]
    readduplication_pdf             // channel: [ val(meta), pdf ]
    readduplication_rscript         // channel: [ val(meta), r   ]

    ch_tin                         // channel: [ val(meta), txt ]

    versions = ch_versions          // channel: [ versions.yml ]
}
