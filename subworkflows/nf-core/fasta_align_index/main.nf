#!/usr/bin/env nextflow

//
// FASTA_ALIGN_INDEX: Align fastq files to a reference genome
//

include { BOWTIE2_BUILD                     } from '../../../modules/nf-core/bowtie2/build/main'
include { BWA_INDEX as BWAMEM1_INDEX        } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                     } from '../../../modules/nf-core/bwamem2/index/main'
include { DRAGMAP_HASHTABLE                 } from '../../../modules/nf-core/dragmap/hashtable/main'
include { SNAPALIGNER_INDEX as SNAP_INDEX   } from '../../../modules/nf-core/snapaligner/index/main'

workflow FASTA_ALIGN_INDEX {

    take:
        ch_fasta        // channel: [ val(meta), fasta ]
        ch_altliftover  // channel: [ val(meta), altliftover ]
        aligner         // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]

    main:

    ch_aligner_index    = Channel.empty()
    ch_versions         = Channel.empty()

    switch (aligner) {
        case 'bowtie2':
            BOWTIE2_BUILD(ch_fasta )                                             // if aligner is bowtie2
            ch_aligner_index = ch_aligner_index.mix(BOWTIE2_BUILD.out.index)
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
            break
        case 'bwamem':
            BWAMEM1_INDEX(ch_fasta)                                              // If aligner is bwa-mem
            // TODO: add altliftover file to index files
            ch_aligner_index = ch_aligner_index.mix(BWAMEM1_INDEX.out.index)
            ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
            break
        case 'bwamem2':
            BWAMEM2_INDEX(ch_fasta)                                              // If aligner is bwa-mem2
            // TODO: add altliftover file to index files
            ch_aligner_index = ch_aligner_index.mix(BWAMEM2_INDEX.out.index)
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
            break
        case 'dragmap':
            DRAGMAP_HASHTABLE(ch_fasta)                                          // If aligner is dragmap
            ch_aligner_index = ch_aligner_index.mix(DRAGMAP_HASHTABLE.out.hashmap)
            ch_versions = ch_versions.mix(DRAGMAP_HASHTABLE.out.versions)
            break
        case 'snap':
            ch_snap_reference = ch_fasta.join(ch_altliftover)
                                .map {meta, fasta, alt -> [meta, fasta, [], [], alt}

            SNAP_INDEX(ch_snap_reference)                                        // If aligner is snap
            ch_aligner_index = ch_aligner_index.mix(SNAP_INDEX.out.index)
            ch_versions = ch_versions.mix(SNAP_INDEX.out.versions)
            break
        default:
            exit 1, "Unknown aligner: ${aligner}"
    }
    ch_aligner_index.dump(tag: 'Aligner Index files')

    emit:
        index    = ch_aligner_index // channel: [ val(meta), index ]
        versions = ch_versions      // channel: [ versions.yml ]
}
