#!/usr/bin/env nextflow

//
// FASTQ_CONTAM_SEQTK_KRAKEN: Subsample FASTQs and perform contamination screening
//
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from "$moduleDir/modules/nf-core/kraken2/kraken2/main"
include { SEQTK_SAMPLE               } from "$moduleDir/modules/nf-core/seqtk/sample/main"

workflow FASTQ_CONTAM_SEQTK_KRAKEN {

    take:
        ch_reads    //channel: [mandatory] meta,reads
        sample_size //string:  [mandatory] number of reads to subsample
        kraken2_db  //string:  [mandatory] path to Kraken2 DB to use for screening

    main:
        ch_reports  = Channel.empty()
        ch_versions = Channel.empty()

        SEQTK_SAMPLE(ch_reads, sample_size)

        ch_versions.mix(SEQTK_SAMPLE.out.versions)

        KRAKEN2(SEQTK_SAMPLE.out.reads,
                kraken2_db,
                false,
                false
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
        ch_reports.mix(KRAKEN2.out.report)

    emit:
        reports  = ch_reports     // channel: [ [meta], log  ]
        versions = ch_versions    // channel: [ versions.yml ]
}
