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

    // Handle different aligners using conditional logic
    if (val_aligner == 'bowtie2') {
        BOWTIE2_BUILD(ch_fasta)
        ch_aligner_index = ch_aligner_index.mix(BOWTIE2_BUILD.out.index)
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    } else if (val_aligner == 'bwamem') {
        BWAMEM1_INDEX(ch_fasta)

        // using index_ inside the map else nf-test collects full file paths
        ch_aligner_index = ch_aligner_index
            .mix(
                BWAMEM1_INDEX.out.index
                .join(ch_altliftover)
                .map { meta, index_, alt -> [meta, index_ + alt] }
            )
        ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
    } else if (val_aligner == 'bwamem2') {
        BWAMEM2_INDEX(ch_fasta)
        // using index_ inside the map else nf-test collects full file paths
        ch_aligner_index = ch_aligner_index
            .mix(
                BWAMEM2_INDEX.out.index
                .join(ch_altliftover)
                .map { meta, index_, alt -> [meta, index_ + alt] }
            )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    } else if (val_aligner == 'dragmap') {
        DRAGMAP_HASHTABLE(ch_fasta)
        ch_aligner_index = ch_aligner_index.mix(DRAGMAP_HASHTABLE.out.hashmap)
        ch_versions = ch_versions.mix(DRAGMAP_HASHTABLE.out.versions)
    } else if (val_aligner == 'snap') {
        ch_snap_reference = ch_fasta
            .join(ch_altliftover)
            .map { meta, fasta_, alt -> [meta, fasta_, [], [], alt] }

        SNAP_INDEX(ch_snap_reference)
        ch_aligner_index = SNAP_INDEX.out.index
        ch_versions = SNAP_INDEX.out.versions
    } else {
        error "Unknown aligner: ${val_aligner}"
    }

    emit:
        index    = ch_aligner_index // channel: [ val(meta), path(index) ]
        versions = ch_versions      // channel: [ path(versions.yml) ]
}
