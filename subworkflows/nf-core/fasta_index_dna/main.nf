#!/usr/bin/env nextflow

//
// FASTA_INDEX_DNA: Build aligner specific index for fasta files
//

include { BOWTIE2_BUILD                     } from '../../../modules/nf-core/bowtie2/build/main'
include { BWA_INDEX as BWAMEM1_INDEX        } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                     } from '../../../modules/nf-core/bwamem2/index/main'
include { DRAGMAP_HASHTABLE                 } from '../../../modules/nf-core/dragmap/hashtable/main'
include { SNAPALIGNER_INDEX as SNAP_INDEX   } from '../../../modules/nf-core/snapaligner/index/main'

workflow FASTA_INDEX_DNA {

    take:
        ch_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_altliftover  // channel: [mandatory, if aligner is bwamem or bwamem2 or snap] [ val(meta), path(altliftover) ]
        val_aligner     // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]

    main:

    ch_aligner_index    = Channel.empty()
    ch_versions         = Channel.empty()

    switch (val_aligner) {
        case 'bowtie2':
            BOWTIE2_BUILD(ch_fasta )                                             // if aligner is bowtie2
            ch_aligner_index = ch_aligner_index.mix(BOWTIE2_BUILD.out.index)
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
            break
        case 'bwamem':
            BWAMEM1_INDEX(ch_fasta)                                              // If aligner is bwa-mem
            ch_aligner_index = ch_aligner_index
                .mix(
                    BWAMEM1_INDEX.out.index
                    .join(ch_altliftover)
                    .map{meta, index, alt -> [meta, index + alt]}
                )

            ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
            break
        case 'bwamem2':
            BWAMEM2_INDEX(ch_fasta)                                              // If aligner is bwa-mem2
            ch_aligner_index = ch_aligner_index
                .mix(
                    BWAMEM2_INDEX.out.index
                    .join(ch_altliftover)
                    .map{meta, index, alt -> [meta, index + alt]}
                )

            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
            break
        case 'dragmap':
            DRAGMAP_HASHTABLE(ch_fasta)                                          // If aligner is dragmap
            ch_aligner_index = ch_aligner_index.mix(DRAGMAP_HASHTABLE.out.hashmap)
            ch_versions = ch_versions.mix(DRAGMAP_HASHTABLE.out.versions)
            break
        case 'snap':
            ch_snap_reference = ch_fasta
                .join(ch_altliftover)
                .map {meta, fasta, alt -> [meta, fasta, [], [], alt]}

            SNAP_INDEX(ch_snap_reference)                                        // If aligner is snap
            ch_aligner_index = ch_aligner_index.mix(SNAP_INDEX.out.index)
            ch_versions = ch_versions.mix(SNAP_INDEX.out.versions)
            break
        default:
            error "Unknown aligner: ${val_aligner}"
    }

    emit:
        index    = ch_aligner_index // channel: [ val(meta), path(index) ]
        versions = ch_versions      // channel: [ path(versions.yml) ]
}
