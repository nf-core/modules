#!/usr/bin/env nextflow

//
// FASTQ_CONTAM_SEQTK_KRAKEN: Subsample FASTQs and perform contamination screening
//
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../../modules/nf-core/kraken2/kraken2/main'
include { SEQTK_SAMPLE               } from '../../../modules/nf-core/seqtk/sample/main'

workflow FASTQ_CONTAM_SEQTK_KRAKEN {

    take:
        ch_reads    //channel: [mandatory] meta,reads
        sample_size //string:  [mandatory] number of reads to subsample
        kraken2_db  //string:  [mandatory] path to Kraken2 DB to use for screening

    main:
        ch_reports  = Channel.empty()
        ch_versions = Channel.empty()

        // Combine all combinations of reads with sample_size(s).
        // Note using more than 1 sample_size can cause file collisions
        // We add n_reads to meta to avoid collisions
        ch_reads
            .combine(sample_size)
            .map{ it ->
                def meta2 = it[0] + [n_reads: it[2]]
                [ meta2, it[1], it[2] ]
            }
            .set { ch_reads_with_n }

        SEQTK_SAMPLE(ch_reads_with_n)

        ch_versions.mix(SEQTK_SAMPLE.out.versions)

        KRAKEN2(SEQTK_SAMPLE.out.reads,
                kraken2_db,
                false,
                false
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
        ch_reports  = ch_reports.mix(KRAKEN2.out.report)

    emit:
        reports  = ch_reports     // channel: [ [meta], log  ]
        versions = ch_versions    // channel: [ versions.yml ]
}
