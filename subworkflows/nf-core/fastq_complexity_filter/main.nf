include { PRINSEQPLUSPLUS } from '../../../modules/nf-core/prinseqplusplus/main'
include { BBMAP_BBDUK     } from '../../../modules/nf-core/bbmap/bbduk/main'
include { FASTP           } from '../../../modules/nf-core/fastp/main'

workflow FASTQ_COMPLEXITY_FILTER {

    take:
    ch_reads                   // channel: [ val(meta), [ fastq ] ]
    val_complexity_filter_tool // string:  [mandatory] tool_name // choose from: ["prinseqplusplus", "bbduk", "fastp"]

    main:

    ch_log           = channel.empty()
    ch_report        = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_versions      = channel.empty()

    if (val_complexity_filter_tool == "prinseqplusplus") {
        PRINSEQPLUSPLUS( ch_reads )

        ch_filtered_reads = PRINSEQPLUSPLUS.out.good_reads
        ch_log            = PRINSEQPLUSPLUS.out.log
        ch_versions       = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())
    } else if (val_complexity_filter_tool == "bbduk") {
        BBMAP_BBDUK( ch_reads, [] )

        ch_filtered_reads = BBMAP_BBDUK.out.reads
        ch_multiqc_files  = ch_multiqc_files.mix(BBMAP_BBDUK.out.log)
        ch_versions       = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
    } else if (val_complexity_filter_tool == "fastp") {
        FASTP(
            ch_reads.map { meta, files ->  [ meta, files, [] ] },
            false,
            false,
            false
        )

        ch_filtered_reads = FASTP.out.reads
        ch_log            = FASTP.out.log
        ch_report         = FASTP.out.html
        ch_multiqc_files  = ch_multiqc_files.mix(FASTP.out.json)
    } else {
        error('Please choose one of the available complexity filtering tools: ["prinseqplusplus", "bbduk", "fastp"]')
    }

    emit:
    filtered_reads = ch_filtered_reads // channel: [ val(meta), [ fastq.gz ] ]
    logfile        = ch_log            // channel: [ val(meta), [ {txt} ] ]
    report         = ch_report         // channel: [ val(meta), [ {html} ] ]
    multiqc_files  = ch_multiqc_files
    versions       = ch_versions       // channel: [ versions.yml ]
}
