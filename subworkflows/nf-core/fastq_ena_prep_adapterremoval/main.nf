include { ADAPTERREMOVAL as ENA_PREP_ADAPTERREMOVAL } from '../../../modules/nf-core/adapterremoval/main'
include { MD5SUM                                    } from '../../../modules/nf-core/md5sum/main'
include { CAT_CAT                                   } from '../../../modules/nf-core/cat/cat/main'

workflow FASTQ_ENA_PREP_ADAPTERREMOVAL {

    take:
    ch_reads        // channel: [ val(meta), [ fastq ], val(phred_plus) ]
    ch_adapter_list // channel [ path(adapter_list) ]

    main:
    ch_fastqs_for_ena = Channel.empty()
    ch_checksum_file  = Channel.empty()
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    // Move phred value to meta map
    ch_ar_input = ch_reads.map {
        meta, fastq, phred_plus ->
            meta2 = meta + [ phred: phred_plus ]
        [ meta2, fastq ]
    }

    ENA_PREP_ADAPTERREMOVAL( ch_ar_input, ch_adapter_list )
    ch_versions      = ch_versions.mix(ENA_PREP_ADAPTERREMOVAL.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(ENA_PREP_ADAPTERREMOVAL.out.settings)

    // Keep appropriate adapter-trimmed output for each sequencing type
    ch_pe_trimmed_fastqs = ENA_PREP_ADAPTERREMOVAL.out.paired_truncated
        .filter {
            meta, fastqs ->
                meta.single_end == false
        }

    ch_pe_trimmed_fastqs_split = ch_pe_trimmed_fastqs
        .multiMap {
            meta, fastqs ->
            r1: [ meta, fastqs[0] ]
            r2: [ meta, fastqs[1] ]
        }

    ch_se_trimmed_fastqs = ENA_PREP_ADAPTERREMOVAL.out.singles_truncated
        .filter {
            meta, fastqs ->
                meta.single_end == true
        }

    ch_fastqs_for_ena = ch_se_trimmed_fastqs.mix(ch_pe_trimmed_fastqs)

    // MD5SUM accepts one file per job.
    ch_fastqs_for_checksum = ch_se_trimmed_fastqs
        .mix(ch_pe_trimmed_fastqs_split.r1)
        .mix(ch_pe_trimmed_fastqs_split.r2)

    MD5SUM ( ch_fastqs_for_checksum )
    ch_versions       = ch_versions.mix(MD5SUM.out.versions.first())

    ch_checksums_to_cat = MD5SUM.out.checksum
        .map{ [ [], it[1] ] }
        .groupTuple( by:0 )

    CAT_CAT ( ch_checksums_to_cat )
    ch_versions       = ch_versions.mix(CAT_CAT.out.versions.first())
    ch_checksum_file  = CAT_CAT.out.file_out


    emit:
    fastq    = ch_fastqs_for_ena // channel: [ val(meta), [ fastq ] ]
    checksum = ch_checksum_file  // channel: [ md5 ]
    multiqc  = ch_multiqc_files  // channel: [ val(meta), path(settings) ]
    versions = ch_versions       // channel: [ versions.yml ]
}

