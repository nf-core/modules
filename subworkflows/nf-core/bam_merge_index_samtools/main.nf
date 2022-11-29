#!/usr/bin/env nextflow

include { SAMTOOLS_INDEX as INDEX } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE } from '../../../modules/nf-core/samtools/merge/main'

workflow BAM_MERGE_INDEX_SAMTOOLS {
    take:
        ch_bam        // channel: [mandatory] [meta, [ bam, ...]
        fasta         // channel: [mandatory] fasta
        fasta_fai     // channel: [mandatory] fasta_fai

    main:
    ch_versions = Channel.empty()

    ch_bam_to_merge = ch_bam.branch{ meta, bam ->
        single:   (bam instanceof List && bam.size() <= 1 ) || (bam instanceof String)
        multiple: bam instanceof List && bam.size() > 1
        other:    true
    }

    MERGE(ch_bam_to_merge.multiple, fasta, fasta_fai)
    ch_versions = ch_versions.mix(MERGE.out.versions.first())

    ch_bam_to_index = ch_bam_to_merge.single
        .mix(MERGE.out.bam)
        .mix(MERGE.out.cram)
    INDEX(ch_bam_to_merge.single.mix(MERGE.out.bam))
    ch_versions = ch_versions.mix(INDEX.out.versions.first())

    ch_bam_index = ch_bam_to_merge.single.map{meta, bam -> [meta, bam[0]]}
        .mix(MERGE.out.bam)
        .join(INDEX.out.index)

    emit:
        bam_index = ch_bam_index
        versions  = ch_versions
}

