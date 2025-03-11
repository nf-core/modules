#!/usr/bin/env nextflow

//
// FASTQ_ALIGN_DNA: Align fastq files to a reference genome
//


include { BOWTIE2_ALIGN                     } from "../../../modules/nf-core/bowtie2/align/main"
include { BWA_MEM as BWAMEM1_MEM            } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM as BWAMEM2_MEM        } from '../../../modules/nf-core/bwamem2/mem/main'
include { DRAGMAP_ALIGN                     } from "../../../modules/nf-core/dragmap/align/main"
include { SNAPALIGNER_ALIGN as SNAP_ALIGN   } from '../../../modules/nf-core/snapaligner/align/main'



workflow FASTQ_ALIGN_DNA {
    take:
        ch_reads            // channel: [mandatory] meta, reads
        ch_aligner_index    // channel: [mandatory] aligner index
        ch_fasta            // channel: [mandatory] fasta file
        aligner             // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]
        sort                // boolean: [mandatory] true -> sort, false -> don't sort

    main:

        ch_bam_index    = Channel.empty()
        ch_bam          = Channel.empty()
        ch_reports      = Channel.empty()
        ch_versions     = Channel.empty()

        // Align fastq files to reference genome and (optionally) sort
        if (aligner == 'bowtie2') {
                BOWTIE2_ALIGN(ch_reads, ch_aligner_index, ch_fasta, false, sort) // if aligner is bowtie2
                ch_bam = ch_bam.mix(BOWTIE2_ALIGN.out.bam)
                ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
        }
        else if (aligner == 'bwamem'){
                BWAMEM1_MEM  (ch_reads, ch_aligner_index, ch_fasta, sort)        // If aligner is bwa-mem
                ch_bam = ch_bam.mix(BWAMEM1_MEM.out.bam)
                ch_bam_index = ch_bam_index.mix(BWAMEM1_MEM.out.csi)
                ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions)
        }
        else if (aligner == 'bwamem2'){
                BWAMEM2_MEM  (ch_reads, ch_aligner_index, ch_fasta, sort)       // If aligner is bwa-mem2
                ch_bam = ch_bam.mix(BWAMEM2_MEM.out.bam)
                ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        }
        else if (aligner == 'dragmap'){
                DRAGMAP_ALIGN(ch_reads, ch_aligner_index, ch_fasta, sort)       // If aligner is dragmap
                ch_bam = ch_bam.mix(DRAGMAP_ALIGN.out.bam)
                ch_reports = ch_reports.mix(DRAGMAP_ALIGN.out.log)
                ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions)
        }
        else if (aligner == 'snap'){
            SNAP_ALIGN   (ch_reads, ch_aligner_index)                       // If aligner is snap
            ch_bam = ch_bam.mix(SNAP_ALIGN.out.bam)
            ch_bam_index.mix(SNAP_ALIGN.out.bai)
            ch_versions = ch_versions.mix(SNAP_ALIGN.out.versions)
        }
        else {
            error "Unknown aligner: ${aligner}"
        }

    emit:
        bam         = ch_bam        // channel: [ [meta], bam       ]
        bam_index   = ch_bam_index  // channel: [ [meta], csi/bai   ]
        reports     = ch_reports    // channel: [ [meta], log       ]
        versions    = ch_versions   // channel: [ versions.yml      ]
}
