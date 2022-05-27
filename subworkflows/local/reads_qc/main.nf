// This subworkflow performs quality control when reads are submitted
// as primary input, and serves as comparison of reads quality between
// original and trimmed reads.


include { FASTQC as FASTQC_UNTRIMMED } from '../../../modules/fastqc/main'
include { FASTQC as FASTQC_TRIMMED   } from '../../../modules/fastqc/main'


workflow READS_QC {

    take:
    untrimmed_reads  // channel: [mandatory] [ val(meta), [ reads ] ]
    trimmed_reads    // channel: [mandatory] [ val(meta), [ reads ] ]


    main:
    ch_versions = Channel.empty()

    FASTQC_UNTRIMMED (
        untrimmed_reads
    )
    ch_versions = ch_versions.mix(FASTQC_UNTRIMMED.out.versions.first())


    FASTQC_TRIMMED (
        trimmed_reads
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())


    emit:
    versions             = ch_versions                  // channel: [ versions.yml ]
    fastqc_untrimmed     = FASTQC_UNTRIMMED.out.zip     // channel: [ val(meta), [ zip ] ]
    fastqc_trimmed       = FASTQC_TRIMMED.out.zip       // channel: [ val(meta), [ zip ] ]


}

