// ok for single end adapter removal
include { TRIMMOMATIC                         } from '../../../modules/nf-core/trimmomatic/main' // both SE and PE
include { CUTADAPT                            } from '../../../modules/nf-core/cutadapt/main'    // both SE and PE
include { TRIMGALORE                          } from '../../../modules/nf-core/trimgalore/main'  // both SE and PE
include { BBMAP_BBDUK                         } from '../../../modules/nf-core/bbmap/bbduk/main' // both SE and PE
// // allows merging of paired end reads, but will work for single end reads as well
include { FASTP                               } from '../../../modules/nf-core/fastp/main'       // both SE and PE + merge
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SE } from '../../../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PE } from '../../../modules/nf-core/adapterremoval/main'
// include { LEEHOM        } from '../../../modules/nf-core/leehom/main'
// // requires paired end because of merging
// include { NGMERGE       } from '../../../modules/nf-core/ngmerge/main'

workflow FASTQ_REMOVEADAPTERS_MERGE {

    take:
    reads                       // channel: [ val(meta), [ reads ] ]
    skip_trimmomatic            // boolean
    skip_cutadapt               // boolean
    skip_trimgalore             // boolean
    skip_bbduk                  // boolean
    contaminants                // channel: [ reads ] // fasta, adapters to remove
    skip_fastp                  // boolean
    fastp_discard_trimmed_pass  // boolean
    fastp_save_trimmed_fail     // boolean
    save_merged                 // boolean
    skip_adapterremoval         // boolean
    text_adapters               // channel: [ txt ] // adapters to remove, in adapterremoval text format
    // skip_leehom              // boolean
    // skip_ngmerge             // boolean

    main:

    ch_reads                             = reads
    ch_merged_reads                      = channel.empty()
    ch_trimmomatic_unpaired_reads        = channel.empty()
    ch_trimmomatic_trim_log              = channel.empty()
    ch_trimmomatic_summary               = channel.empty()
    ch_trimgalore_log                    = channel.empty()
    ch_trimgalore_unpaired               = channel.empty()
    ch_trimgalore_html                   = channel.empty()
    ch_trimgalore_zip                    = channel.empty()
    ch_fastp_html                        = channel.empty()
    ch_fastp_log                         = channel.empty()
    ch_fastp_reads_fail                  = channel.empty()
    ch_adapterremoval_discarded          = channel.empty()
    ch_adapterremoval_paired_interleaved = channel.empty()
    ch_versions                          = channel.empty()
    ch_multiqc_files                     = channel.empty()

    if (!skip_trimmomatic) {
        TRIMMOMATIC( ch_reads )
        ch_reads                      = TRIMMOMATIC.out.trimmed_reads
        ch_trimmomatic_unpaired_reads = TRIMMOMATIC.out.unpaired_reads
        ch_trimmomatic_trim_log       = TRIMMOMATIC.out.trim_log
        ch_trimmomatic_summary        = TRIMMOMATIC.out.summary
        ch_versions                   = ch_versions.mix(TRIMMOMATIC.out.versions.first())
        ch_multiqc_files              = ch_multiqc_files.mix(TRIMMOMATIC.out.out_log)
    }

    if (!skip_cutadapt) {
        CUTADAPT( ch_reads )
        ch_reads         = CUTADAPT.out.reads
        ch_versions      = ch_versions.mix(CUTADAPT.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log)
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
        ch_reads         = BBMAP_BBDUK.out.reads
        ch_versions      = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log)
    }

    if (!skip_fastp) {
        FASTP(
            ch_reads.map { meta, files ->  [ meta, files, contaminants ] },
            fastp_discard_trimmed_pass,
            fastp_save_trimmed_fail,
            save_merged
        )
        ch_reads              = FASTP.out.reads
        ch_merged_reads       = FASTP.out.reads_merged
        ch_fastp_html         = FASTP.out.html
        ch_fastp_log          = FASTP.out.log
        ch_fastp_reads_fail   = FASTP.out.reads_fail
        ch_versions           = ch_versions.mix(FASTP.out.versions.first())
        ch_multiqc_files      = ch_multiqc_files.mix(FASTP.out.json)
    }

    if (!skip_adapterremoval) {
        ch_adapterremoval_in = ch_reads
            .branch { meta, _reads ->
                single: meta.single_end
                paired: !meta.single_end
            }

        ADAPTERREMOVAL_SE(ch_adapterremoval_in.single, text_adapters)
        ch_versions = ch_versions.mix(ADAPTERREMOVAL_SE.out.versions.first())
        ADAPTERREMOVAL_PE(ch_adapterremoval_in.paired, text_adapters)
        ch_versions = ch_versions.mix(ADAPTERREMOVAL_PE.out.versions.first())

        ch_reads        = ADAPTERREMOVAL_SE.out.singles_truncated.mix(ADAPTERREMOVAL_PE.out.paired_truncated)
        ch_merged_reads = ADAPTERREMOVAL_SE.out.collapsed
            .mix(
                ADAPTERREMOVAL_PE.out.collapsed,
                ADAPTERREMOVAL_SE.out.collapsed_truncated,
                ADAPTERREMOVAL_PE.out.collapsed_truncated
            )
        ch_adapterremoval_discarded           = ADAPTERREMOVAL_SE.out.discarded.mix(ADAPTERREMOVAL_PE.out.discarded)
        ch_adapterremoval_paired_interleaved  = ADAPTERREMOVAL_SE.out.paired_interleaved.mix(ADAPTERREMOVAL_PE.out.paired_interleaved)
        ch_multiqc_files                      = ch_multiqc_files.mix(ADAPTERREMOVAL_PE.out.settings, ADAPTERREMOVAL_SE.out.settings)
    }

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
    reads                              = ch_reads                               // channel: [ val(meta), [ fastq.gz ] ]
    merged_reads                       = ch_merged_reads                        // channel: [ val(meta), [ fastq.gz ] ]
    trimmomatic_unpaired_reads         = ch_trimmomatic_unpaired_reads          // channel: [ val(meta), [ fastq.gz ] ]
    trimmomatic_trim_log               = ch_trimmomatic_trim_log                // channel: [ val(meta), [ log ] ]
    trimmomatic_summary                = ch_trimmomatic_summary                 // channel: [ val(meta), [ summary ] ]
    trimgalore_log                     = ch_trimgalore_log                      // channel: [ val(meta), [ txt ] ]
    trimgalore_unpaired                = ch_trimgalore_unpaired                 // channel: [ val(meta), [ fq.gz ] ]
    trimgalore_html                    = ch_trimgalore_html                     // channel: [ val(meta), [ html ] ]
    trimgalore_zip                     = ch_trimgalore_zip                      // channel: [ val(meta), [ zip ] ]
    fastp_html                         = ch_fastp_html                          // channel: [ val(meta), [ html ] ]
    fastp_log                          = ch_fastp_log                           // channel: [ val(meta), [ log ] ]
    fastp_reads_fail                   = ch_fastp_reads_fail                    // channel: [ val(meta), [ fastq.gz ] ]
    adapterremoval_discarded           = ch_adapterremoval_discarded            // channel: [ val(meta), [ fastq.gz ] ]
    adapterremoval_paired_interleaved  = ch_adapterremoval_paired_interleaved   // channel: [ val(meta), [ fastq.gz ] ]
    versions                           = ch_versions                            // channel: [ versions.yml ]
    multiqc_files                      = ch_multiqc_files
}
