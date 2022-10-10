#!/usr/bin/env nextflow

//
// FASTQ_ALIGN: Align fastq files to a reference genome
//


include { BOWTIE2_ALIGN                     } from "../../../modules/nf-core/bowtie2/align/main"
include { BWA_MEM as BWAMEM1_MEM            } from '../../../modules/nf-core/bwa/mem/main'
include { BWA_MEM2 as BWAMEM2_MEM           } from '../../../modules/nf-core/bwamem2/mem/main'
include { DRAGMAP_ALIGN                     } from "../../../modules/nf-core/dragmap/align/main"
include { SNAPALIGNER_ALIGN as SNAP_ALIGN   } from '../../../modules/nf-core/snapaligner/align/main'



workflow FASTQ_ALIGN {
    take:
        ch_reads            // channel: [mandatory] meta, reads
        ch_aligner_index    // channel: [mandatory] aligner index
        aligner             // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]
        sort                // boolean: [mandatory] true -> sort, false -> don't sort

    main:

        ch_versions = Channel.empty()

        // Align fastq files to reference genome and (optionally) sort
        switch (aligner) {
            case 'bowtie2':
                BOWTIE2_ALIGN(ch_reads, ch_aligner_index, false, sort) // if aligner is bowtie2
                break
            case 'bwa_mem':
                BWAMEM1_MEM  (ch_reads, ch_aligner_index, sort)        // If aligner is bwa-mem
                break
            case 'bwa_mem2':
                BWAMEM2_MEM  (ch_reads, ch_aligner_index, sort)        // If aligner is bwa-mem2
                break
            case 'dragmap':
                DRAGMAP_ALIGN(ch_reads, ch_map_index, sort)            // If aligner is dragmap
                break
            case 'snap':
                SNAP_ALIGN   (ch_reads, ch_aligner_index)              // If aligner is snap
                break
            default:
                exit 1, "Unknown aligner: ${aligner}"
        }

        ch_versions = ch_versions.mix(
            BOWTIE2_ALIGN.out.versions,
            BWAMEM1_MEM.out.versions,
            BWA_MEM2_MEM.out.versions,
            DRAGMAP_ALIGN.out.versions,
            SNAP_ALIGN.out.versions
        )

        // Get the bam files from the aligner
        // Only one aligner is run
        ch_bam = Channel.empty()
        ch_bam = ch_bam.mix(
            BOWTIE2_ALIGN.out.bam,
            BWAMEM1_MEM.out.bam,
            BWA_MEM2_MEM.out.bam,
            DRAGMAP_ALIGN.out.bam,
            SNAP_ALIGN.out.bam,
        )
        ch_cleaned_bam = ch_bam.map { meta, bam ->
            [ meta.findAll { !(it.key in ['readgroup', 'library']) }, bam ]
        }.groupTuple()

        // Gather indexes where available
        ch_bai = Channel.empty()
        ch_bai= ch_bai.mix(
            SNAP_ALIGN.out.bai,
        )
        ch_cleaned_bai = ch_bai.map { meta, bai ->
            [ meta.findAll { !(it.key in ['readgroup', 'library']) }, bai ]
        }.groupTuple()

        // Gather reports where available
        ch_reports = ch_reports.mix(DRAGMAP_ALIGN.out.log)


    emit:
        bam      = ch_cleaned_bam // channel: [ [meta], bam  ]
        bai      = ch_cleaned_bai // channel: [ [meta], bai  ]
        reports  = ch_reports     // channel: [ [meta], log  ]
        versions = ch_versions    // channel: [ versions.yml ]
}
