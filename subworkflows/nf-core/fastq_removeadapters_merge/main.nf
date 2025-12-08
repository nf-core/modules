// ok for single end adapter removal
include { TRIMMOMATIC   } from '../../../modules/nf-core/trimmomatic/main' // both SE and PE
include { CUTADAPT      } from '../../../modules/nf-core/cutadapt/main'    // both SE and PE
include { TRIMGALORE    } from '../../../modules/nf-core/trimgalore/main'  // both SE and PE
include { BBMAP_BBDUK   } from '../../../modules/nf-core/bbmap/bbduk/main' // both SE and PE
// // allows merging of paired end reads, but will work for single end reads as well
// include { FASTP         } from '../../../modules/nf-core/fastp/main'
// include { ADAPTERREMOVAL} from '../../../modules/nf-core/adapterremoval/main'
// include { LEEHOM        } from '../../../modules/nf-core/leehom/main'
// // requires paired end because of merging
// include { NGMERGE       } from '../../../modules/nf-core/ngmerge/main'

workflow FASTQ_REMOVEADAPTERS_MERGE {

    take:
    reads                // channel: [ val(meta), [ reads ] ]
    skip_trimmomatic     // boolean
    skip_cutadapt        // boolean
    skip_trimgalore      // boolean
    skip_bbduk           // boolean
    contaminants         // channel: [ reads ]
    // skip_fastp           // boolean
    // skip_adapterremoval  // boolean
    // skip_leehom          // boolean
    // skip_ngmerge         // boolean
    // do_merge
    // adapters_

    main:

    ch_reads                      = reads
    ch_trimmomatic_unpaired_reads = channel.empty()
    ch_trimmomatic_trim_log       = channel.empty()
    ch_trimmomatic_out_log        = channel.empty()
    ch_trimmomatic_summary        = channel.empty()
    ch_cutadapt_log               = channel.empty()
    ch_trimgalore_log             = channel.empty()
    ch_trimgalore_unpaired        = channel.empty()
    ch_trimgalore_html            = channel.empty()
    ch_trimgalore_zip             = channel.empty()
    ch_bbduk_log                  = channel.empty()
    ch_versions                   = channel.empty()

    if (!skip_trimmomatic) {
        TRIMMOMATIC( ch_reads )
        ch_reads                      = TRIMMOMATIC.out.trimmed_reads
        ch_trimmomatic_unpaired_reads = TRIMMOMATIC.out.unpaired_reads
        ch_trimmomatic_trim_log       = TRIMMOMATIC.out.trim_log
        ch_trimmomatic_out_log        = TRIMMOMATIC.out.out_log
        ch_trimmomatic_summary        = TRIMMOMATIC.out.summary
        ch_versions                   = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    }

    if (!skip_cutadapt) {
        CUTADAPT( ch_reads )
        ch_reads        = CUTADAPT.out.reads
        ch_cutadapt_log = CUTADAPT.out.log
        ch_versions     = ch_versions.mix(CUTADAPT.out.versions.first())
    }

    if (!skip_trimgalore) {
        TRIMGALORE( ch_reads )
        ch_reads               = TRIMGALORE.out.reads
        ch_trimgalore_log      = TRIMGALORE.out.log
        ch_trimgalore_unpaired = TRIMGALORE.out.unpaired
        ch_trimgalore_html     = TRIMGALORE.out.html
        ch_trimgalore_zip      = TRIMGALORE.out.zip
        ch_versions            = ch_versions.mix(TRIMGALORE.out.versions.first())
    }

    if (!skip_bbduk) {
        BBMAP_BBDUK( ch_reads, contaminants )
        ch_reads     = BBMAP_BBDUK.out.reads
        ch_bbduk_log = BBMAP_BBDUK.out.log
        ch_versions  = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
    }

//     if (!skip_fastp && !do_merge) {
//         ch_reads.view { "DEBUG: BEFORE FASTP → $it" }
//         FASTP( ch_reads.map { meta, r ->  tuple(meta, r, adapters_contaminants)},
//             false,
//             true,
//             do_merge)
//         ch_reads = FASTP.out.reads.map { meta, r -> [meta, r] }
//         ch_versions = ch_versions.mix(FASTP.out.versions)
//         ch_reads.view { "DEBUG: AFTER FASTP → $it" }
//     }

//     if (!skip_adapterremoval && !do_merge) {
//     ADAPTERREMOVAL( ch_reads, adapters_contaminants )
//     ch_reads.view { "DEBUG: BEFORE ADAPTERREMOVAL → $it" }
//          ch_reads = ADAPTERREMOVAL.out.paired_truncated
//     .mix(ADAPTERREMOVAL.out.singles_truncated)
//     .map { tup ->
//         def meta  = tup[0]
//         def r = tup[1]
//         [meta, r]
//     }
//          ch_versions = ch_versions.mix(ADAPTERREMOVAL.out.versions)
//          ch_reads.view { "DEBUG: AFTER ADAPTERREMOVAL → $it" }
//     }

//     if (!skip_leehom && !do_merge) {
//         ch_reads.view { "DEBUG: BEFORE LEEHOM → $it" }
//         ch_leehom_input = ch_reads.map { meta, r ->
//         if (meta.single_end)
//             return [meta, r]
//         else
//             return [meta, [reads[0], reads[1]]]
//     }

//     LEEHOM(ch_leehom_input)

//         ch_leehom_pe =
//     LEEHOM.out.unmerged_r1_fq_pass
//         .combine(LEEHOM.out.unmerged_r2_fq_pass)
//         .map { left, right ->
//             def meta = left[0]
//             def r1   = left[1]
//             def r2   = right[1]
//             [meta, [r1, r2]]
//         }
//         .filter { meta, r -> !meta.single_end }  // only PE

// // ---- SE OUTPUT ----
// ch_leehom_se =
//     LEEHOM.out.fq_pass
//         .filter { meta, r -> meta.single_end }
//         .map { meta, r -> [meta, r] }

// // ---- MERGE CHANNELS ----
// ch_reads = ch_leehom_pe.mix(ch_leehom_se)

//         ch_versions = ch_versions.mix(LEEHOM.out.versions)
//     }

    emit:
    reads                      = ch_reads                       // channel: [ val(meta), [ fastq.gz ] ]
    trimmomatic_unpaired_reads = ch_trimmomatic_unpaired_reads  // channel: [ val(meta), [ fastq.gz ] ]
    trimmomatic_trim_log       = ch_trimmomatic_trim_log        // channel: [ val(meta), [ log ] ]
    trimmomatic_out_log        = ch_trimmomatic_out_log         // channel: [ val(meta), [ log ] ]
    trimmomatic_summary        = ch_trimmomatic_summary         // channel: [ val(meta), [ summary ] ]
    cutadapt_log               = ch_cutadapt_log                // channel: [ val(meta), [ log ] ]
    trimgalore_log             = ch_trimgalore_log              // channel: [ val(meta), [ txt ] ]
    trimgalore_unpaired        = ch_trimgalore_unpaired         // channel: [ val(meta), [ fq.gz ] ]
    trimgalore_html            = ch_trimgalore_html             // channel: [ val(meta), [ html ] ]
    trimgalore_zip             = ch_trimgalore_zip              // channel: [ val(meta), [ zip ] ]
    bbduk_log                  = ch_bbduk_log                   // channel: [ val(meta), [ log ] ]
    versions                   = ch_versions                    // channel: [ versions.yml ]
}
