// both SE and PE
include { TRIMMOMATIC                         } from '../../../modules/nf-core/trimmomatic/main'
include { CUTADAPT                            } from '../../../modules/nf-core/cutadapt/main'
include { TRIMGALORE                          } from '../../../modules/nf-core/trimgalore/main'
include { BBMAP_BBDUK                         } from '../../../modules/nf-core/bbmap/bbduk/main'
include { LEEHOM                              } from '../../../modules/nf-core/leehom/main'
// both SE and PE, plus merging
include { FASTP                               } from '../../../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SE } from '../../../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PE } from '../../../modules/nf-core/adapterremoval/main'
// helper module for concatenating adapterremoval paired-end processed reads
include { CAT_FASTQ                           } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_REMOVEADAPTERS_MERGE {

    take:
    ch_input_reads                 // channel: [mandatory] meta, reads
    val_adapter_tool               // string:  [mandatory] tool_name // choose from: ["trimmomatic", "cutadapt", "trimgalore", "bbduk", "leehom", "fastp", "adapterremoval"]
    ch_custom_adapters_file        // channel: [optional]  {fasta,txt} // fasta, for bbduk or fastp, or txt, for adapterremoval
    val_save_merged                // boolean: [mandatory] if true, will return the merged reads instead, for fastp and adapterremoval
    val_fastp_discard_trimmed_pass // boolean: [mandatory] // only for fastp
    val_fastp_save_trimmed_fail    // boolean: [mandatory] // only for fastp

    main:

    ch_discarded_reads    = channel.empty() // from trimmomatic, trimgalore, leehom, fastp, adapterremoval
    ch_log                = channel.empty() // from trimmomatic, trimgalore, fastp
    ch_report             = channel.empty() // from trimmomatic, trimgalore, fastp
    ch_versions           = channel.empty()
    ch_multiqc_files      = channel.empty() // from trimmomatic, cutadapt, bbduk, leehom, fastp, adapterremoval

    if (val_adapter_tool == "trimmomatic") {
        TRIMMOMATIC( ch_input_reads )

        ch_processed_reads = TRIMMOMATIC.out.trimmed_reads
        ch_discarded_reads = ch_discarded_reads.mix(TRIMMOMATIC.out.unpaired_reads.transpose()) // .transpose() because paired reads will output 2 unpaired files in an array
        ch_log             = TRIMMOMATIC.out.trim_log
        ch_report          = TRIMMOMATIC.out.summary
        ch_multiqc_files   = ch_multiqc_files.mix(TRIMMOMATIC.out.out_log)
    } else if (val_adapter_tool == "cutadapt") {
        CUTADAPT( ch_input_reads )

        ch_processed_reads = CUTADAPT.out.reads
        ch_multiqc_files   = ch_multiqc_files.mix(CUTADAPT.out.log)
    } else if (val_adapter_tool == "trimgalore") {
        TRIMGALORE( ch_input_reads )

        ch_processed_reads = TRIMGALORE.out.reads
        ch_discarded_reads = ch_discarded_reads.mix(TRIMGALORE.out.unpaired)
        ch_log             = TRIMGALORE.out.log
        ch_report          = TRIMGALORE.out.html.mix(TRIMGALORE.out.zip)
    } else if (val_adapter_tool == "bbduk") {
        BBMAP_BBDUK( ch_input_reads, ch_custom_adapters_file )

        ch_processed_reads = BBMAP_BBDUK.out.reads
        ch_versions        = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        ch_multiqc_files   = ch_multiqc_files.mix(BBMAP_BBDUK.out.log)
    } else if (val_adapter_tool == "leehom") {
        LEEHOM( ch_input_reads )

        ch_processed_reads = LEEHOM.out.fq_pass
            .join(LEEHOM.out.unmerged_r1_fq_pass, by: 0, remainder: true)
            .join(LEEHOM.out.unmerged_r2_fq_pass, by: 0, remainder: true)
            .map { meta, single, r1, r2 ->
                if (meta.single_end) {
                    return [meta, single]
                } else {
                    return [meta, [r1, r2]]
                }
            }
        ch_discarded_reads = ch_discarded_reads.mix(LEEHOM.out.fq_fail, LEEHOM.out.unmerged_r1_fq_fail, LEEHOM.out.unmerged_r2_fq_fail)
        ch_versions        = ch_versions.mix(LEEHOM.out.versions.first())
        ch_multiqc_files   = ch_multiqc_files.mix(LEEHOM.out.log)
    } else if (val_adapter_tool == "fastp") {
        FASTP(
            ch_input_reads.map { meta, files ->  [ meta, files, ch_custom_adapters_file ] },
            val_fastp_discard_trimmed_pass,
            val_fastp_save_trimmed_fail,
            val_save_merged
        )

        if (val_save_merged) {
            ch_processed_reads = FASTP.out.reads_merged
        } else {
            ch_processed_reads = FASTP.out.reads
        }
        ch_discarded_reads = ch_discarded_reads.mix(FASTP.out.reads_fail.transpose()) // .transpose() because paired reads have 3 fail files in an array
        ch_log             = FASTP.out.log
        ch_report          = FASTP.out.html
        ch_multiqc_files   = ch_multiqc_files.mix(FASTP.out.json)
    } else if (val_adapter_tool == "adapterremoval") {
        ch_adapterremoval_in = ch_input_reads
            .branch { meta, _reads ->
                single: meta.single_end
                paired: !meta.single_end
            }

        ADAPTERREMOVAL_SE( ch_adapterremoval_in.single, ch_custom_adapters_file )
        ADAPTERREMOVAL_PE( ch_adapterremoval_in.paired, ch_custom_adapters_file )

        if (val_save_merged) { // merge
            ch_concat_fastq = channel.empty()
                .mix(
                    ADAPTERREMOVAL_PE.out.collapsed,
                    ADAPTERREMOVAL_PE.out.collapsed_truncated,
                    ADAPTERREMOVAL_PE.out.singles_truncated,
                )
                .map { meta, reads ->
                    def meta_new = meta.clone()
                    meta_new.single_end = true
                    [meta_new, reads]
                }
                .groupTuple()
                // Paired-end reads cause a nested tuple during grouping.
                // We want to present a flat list of files to `CAT_FASTQ`.
                .map { meta, fastq -> [meta, fastq.flatten()] }

            CAT_FASTQ( ch_concat_fastq )

            ch_processed_reads = CAT_FASTQ.out.reads.mix(ADAPTERREMOVAL_SE.out.singles_truncated)
        } else { // no merge
            ch_processed_reads = ADAPTERREMOVAL_PE.out.paired_truncated.mix(ADAPTERREMOVAL_SE.out.singles_truncated)
        }
        ch_discarded_reads    = ch_discarded_reads.mix(ADAPTERREMOVAL_SE.out.discarded, ADAPTERREMOVAL_PE.out.discarded)
        ch_versions           = ch_versions.mix(ADAPTERREMOVAL_SE.out.versions.first(), ADAPTERREMOVAL_PE.out.versions.first())
        ch_multiqc_files      = ch_multiqc_files.mix(ADAPTERREMOVAL_PE.out.settings, ADAPTERREMOVAL_SE.out.settings)
    } else {
        error('Please choose one of the available adapter removal and merging tools: ["trimmomatic", "cutadapt", "trimgalore", "bbduk", "leehom", "fastp", "adapterremoval"]')
    }

    emit:
    processed_reads    = ch_processed_reads    // channel: [ val(meta), [ fastq.gz ] ]
    discarded_reads    = ch_discarded_reads    // channel: [ val(meta), [ fastq.gz ] ]
    logfile            = ch_log                // channel: [ val(meta), [ {log,txt} ] ]
    report             = ch_report             // channel: [ val(meta), [ {summary,html,zip} ] ]
    versions           = ch_versions           // channel: [ versions.yml ]
    multiqc_files      = ch_multiqc_files
}
