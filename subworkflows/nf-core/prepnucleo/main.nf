// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FGBIO_FASTQTOBAM } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { PICARD_MERGESAMFILES } from '../../../modules/nf-core/picard/mergesamfiles/'
include { GATK4_SAMTOFASTQ } from '../../../modules/nf-core/gatk4/samtofastq/main'
include { FASTP } from '../../../modules/nf-core/fastp/main'


workflow PREPNUCLEO {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fastq // channel: [ val(meta), [ bam ] ]


    main:

    ch_versions = Channel.empty()

    // FGBIO_FASTQTOBAM: get unmerged bams 
    // ch_fastq is a channel, which enables parallel
    // channels enable parallel: https://www.nextflow.io/docs/latest/faq.html?highlight=parallel
    FGBIO_FASTQTOBAM (
        ch_fastq
    )
    FGBIO_FASTQTOBAM.out.bam.map{
        meta, bam ->
            [bam]
    }.collect().map{
        bams ->
         [[id: 'unmerged_bams'], bams ]
    }.set{unmerged_bams}
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions)

    // PICARD_MERGESAMFILES: merge bams files 
    PICARD_MERGESAMFILES (
        unmerged_bams
    ).bam.map {
        meta, bam ->
            new_id = 'merged_bam'
            [[id: new_id], bam ]
    }.set {merged_bam}
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions)

    // GATK4_SAMTOFASTQ: get fastqs from merged bam 
    GATK4_SAMTOFASTQ (
        merged_bam
    ).fastq.map {
        meta, fastq ->
            new_id = 'merged_fastq'
            [[id: new_id], fastq ]
    }.set {merged_fastq}
    ch_versions = ch_versions.mix(GATK4_SAMTOFASTQ.out.versions)

    // GATK4_SAMTOFASTQ: Run fastp on fastqs 
    FASTP (
        merged_fastq, [], false, false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // final emit
    emit:
    // TODO nf-core: edit emitted channels
    bam      = FASTP.out.reads          

    versions = ch_versions                     // channel: [ versions.yml ]
}


