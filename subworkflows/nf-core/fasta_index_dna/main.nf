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
        ch_altliftover  // channel: [optional, only used if aligner is snap] [ val(meta), path(altliftover) ]
        val_aligner     // string:  [mandatory] aligner [bowtie2, bwamem, bwamem2, dragmap, snap]

    main:

    ch_aligner_index    = channel.empty()

    if (val_aligner != "snap") {
        ch_altliftover
            .filter { _meta, altliftover -> altliftover != [] }
            .subscribe { error "Liftover file is only currently supported for aligner `snap`" }
    }

    // Handle different aligners using conditional logic
    if (val_aligner == 'bowtie2') {
        BOWTIE2_BUILD(ch_fasta)
        ch_aligner_index = BOWTIE2_BUILD.out.index
    } else if (val_aligner == 'bwamem') {
        BWAMEM1_INDEX(ch_fasta)
        ch_aligner_index = BWAMEM1_INDEX.out.index
        // bwa/index uses a topic channel for versions
    } else if (val_aligner == 'bwamem2') {
        BWAMEM2_INDEX(ch_fasta)
        ch_aligner_index = BWAMEM2_INDEX.out.index
    } else if (val_aligner == 'dragmap') {
        DRAGMAP_HASHTABLE(ch_fasta)
        ch_aligner_index = DRAGMAP_HASHTABLE.out.hashmap
    } else if (val_aligner == 'snap') {
        ch_snap_reference = ch_fasta
            .combine(ch_altliftover)
            .map { meta, fasta_, _meta2, liftover -> [meta, fasta_, [], [], liftover] }
        SNAP_INDEX(ch_snap_reference)
        ch_aligner_index = SNAP_INDEX.out.index
    } else {
        error "Unknown aligner: ${val_aligner}"
    }

    emit:
        index    = ch_aligner_index // channel: [ val(meta), path(index) ]
}
