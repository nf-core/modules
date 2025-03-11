include { UNTAR                     } from '../../../modules/nf-core/untar/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX             } from '../../../modules/nf-core/bwameth/index/main'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'

workflow FASTA_INDEX_BISMARK_BWAMETH {

    take:
    fasta         // channel: [ val(meta), [ fasta ] ]
    fasta_index   // channel: [ val(meta), [ fasta index ] ]
    bismark_index // channel: [ val(meta), [ bismark index ] ]
    bwameth_index // channel: [ val(meta), [ bwameth index ] ]

    main:

    ch_fasta         = Channel.empty()
    ch_fasta_index   = Channel.empty()
    ch_bismark_index = Channel.empty()
    ch_bwameth_index = Channel.empty()
    ch_versions      = Channel.empty()

    if (fasta.toString().endsWith('.gz')) {
        GUNZIP (
            [ [:], file(fasta, checkIfExists: true) ]
        )
        ch_fasta    = GUNZIP.out.gunzip
        ch_versions = ch_versions.mix(GUNZIP.out.versions)
    } else {
        ch_fasta    = Channel.value([[:], file(fasta, checkIfExists: true)])
    }

    // Aligner: bismark or bismark_hisat
    if( params.aligner =~ /bismark/ ){
        /*
         * Generate bismark index if not supplied
         */
        if (bismark_index) {
            if (bismark_index.toString().endsWith('.gz')) {
                UNTAR (
                    [ [:], file(bismark_index, checkIfExists: true) ]
                )
                ch_bismark_index = UNTAR.out.untar
                ch_versions      = ch_versions.mix(UNTAR.out.versions)
            } else {
                ch_bismark_index = Channel.value([[:], file(bismark_index, checkIfExists: true)])
            }
        } else {
            BISMARK_GENOMEPREPARATION (
                ch_fasta
            )
            ch_bismark_index = BISMARK_GENOMEPREPARATION.out.index
            ch_versions      = ch_versions.mix(BISMARK_GENOMEPREPARATION.out.versions)
        }
    }

    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' ){
        /*
         * Generate bwameth index if not supplied
         */
        if (bwameth_index) {
            if (bwameth_index.toString().endsWith('.gz')) {
                UNTAR (
                    [ [:], file(bwameth_index, checkIfExists: true) ]
                )
                ch_bwameth_index = UNTAR.out.untar
                ch_versions      = ch_versions.mix(UNTAR.out.versions)
            } else {
                ch_bwameth_index = Channel.value([[:], file(bwameth_index, checkIfExists: true)])
            }
        } else {
            BWAMETH_INDEX (
                ch_fasta
            )
            ch_bwameth_index = BWAMETH_INDEX.out.index
            ch_versions      = ch_versions.mix(BWAMETH_INDEX.out.versions)
        }

        /*
         * Generate fasta index if not supplied
         */
        if (fasta_index) {
            ch_fasta_index = Channel.value(file(fasta_index, checkIfExists: true))
        } else {
            SAMTOOLS_FAIDX(
                ch_fasta,
                [[:], []],
                false // No sizes generation
            )
            ch_fasta_index = SAMTOOLS_FAIDX.out.fai
            ch_versions    = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta         // channel: [ val(meta), [ fasta ] ]
    fasta_index   = ch_fasta_index   // channel: [ val(meta), [ fasta index ] ]
    bismark_index = ch_bismark_index // channel: [ val(meta), [ bismark index ] ]
    bwameth_index = ch_bwameth_index // channel: [ val(meta), [ bwameth index ] ]
    versions      = ch_versions      // channel: [ versions.yml ]
}

