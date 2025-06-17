include { UNTAR                     } from '../../../modules/nf-core/untar/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX             } from '../../../modules/nf-core/bwameth/index/main'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'

workflow FASTA_INDEX_BISMARK_BWAMETH {

    take:
    fasta            // channel: [ val(meta), [ fasta ] ]
    fasta_index      // channel: [ val(meta), [ fasta index ] ]
    bismark_index    // channel: [ val(meta), [ bismark index ] ]
    bwameth_index    // channel: [ val(meta), [ bwameth index ] ]
    aligner          // string: bismark, bismark_hisat or bwameth
    collecthsmetrics // boolean: whether to run picard collecthsmetrics

    main:

    ch_fasta         = Channel.empty()
    ch_fasta_index   = Channel.empty()
    ch_bismark_index = Channel.empty()
    ch_bwameth_index = Channel.empty()
    ch_versions      = Channel.empty()

    // Check if fasta file is gzipped and decompress if needed
    fasta
        .branch {
            gzipped: it[1].toString().endsWith('.gz')
            unzipped: true
        }
        .set { ch_fasta_branched }

    GUNZIP (
        ch_fasta_branched.gzipped.ifEmpty([])
    )

    ch_fasta    = ch_fasta_branched.unzipped.mix(GUNZIP.out.gunzip)
    ch_versions = ch_versions.mix(GUNZIP.out.versions)

    // Aligner: bismark or bismark_hisat
    if( aligner =~ /bismark/ ){
        /*
         * Generate bismark index if not supplied
         */
        if (bismark_index) {
            // Handle channel-based bismark index
            bismark_index
                .branch {
                    gzipped: it[1].toString().endsWith('.gz')
                    unzipped: true
                }
                .set { ch_bismark_index_branched }

            UNTAR (
                ch_bismark_index_branched.gzipped.ifEmpty([])
            )

            ch_bismark_index = ch_bismark_index_branched.unzipped.mix(UNTAR.out.untar)
            ch_versions      = ch_versions.mix(UNTAR.out.versions)
        } else {
            BISMARK_GENOMEPREPARATION (
                ch_fasta
            )
            ch_bismark_index = BISMARK_GENOMEPREPARATION.out.index
            ch_versions      = ch_versions.mix(BISMARK_GENOMEPREPARATION.out.versions)
        }
    }

    // Aligner: bwameth
    else if ( aligner == 'bwameth' ){
        /*
         * Generate bwameth index if not supplied
         */
        if (bwameth_index) {
            // Handle channel-based bwameth index
            bwameth_index
                .branch {
                    gzipped: it[1].toString().endsWith('.gz')
                    unzipped: true
                }
                .set { ch_bwameth_index_branched }

            UNTAR (
                ch_bwameth_index_branched.gzipped.ifEmpty([])
            )

            ch_bwameth_index = ch_bwameth_index_branched.unzipped.mix(UNTAR.out.untar)
            ch_versions      = ch_versions.mix(UNTAR.out.versions)
        } else {
            BWAMETH_INDEX (
                ch_fasta
            )
            ch_bwameth_index = BWAMETH_INDEX.out.index
            ch_versions      = ch_versions.mix(BWAMETH_INDEX.out.versions)
        }
    }

    /*
    * Generate fasta index if not supplied for bwameth workflow or picard collecthsmetrics tool
    */
    if (aligner == 'bwameth' | collecthsmetrics) {
        // already exising fasta index
        if (fasta_index) {
            ch_fasta_index = fasta_index
        } else {
            SAMTOOLS_FAIDX(
                ch_fasta,
                [[:], []],
                false
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

