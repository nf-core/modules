include { UNTAR                     } from '../../../modules/nf-core/untar/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX             } from '../../../modules/nf-core/bwameth/index/main'
include { BWA_INDEX                 } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'

workflow FASTA_INDEX_BISMARK_BWAMETH_BWAMEM {

    take:
    fasta            // channel: [ val(meta), [ fasta ] ]
    fasta_index      // channel: [ val(meta), [ fasta index ] ]
    bismark_index    // channel: [ val(meta), [ bismark index ] ]
    bwameth_index    // channel: [ val(meta), [ bwameth index ] ]
    bwamem_index     // channel: [ val(meta), [ bwamem index ] ]
    aligner          // string: bismark, bismark_hisat or bwameth
    collecthsmetrics // boolean: whether to run picard collecthsmetrics
    use_mem2         // boolean: generate mem2 index if no index provided, and bwameth is selected

    main:

    ch_fasta         = Channel.empty()
    ch_fasta_index   = Channel.empty()
    ch_bismark_index = Channel.empty()
    ch_bwameth_index = Channel.empty()
    ch_bwamem_index  = Channel.empty()
    ch_versions      = Channel.empty()

    // Check if fasta file is gzipped and decompress if needed
    fasta
        .branch {
            gzipped: it[1].toString().endsWith('.gz')
            unzipped: true
        }
        .set { ch_fasta_branched }

    GUNZIP (
        ch_fasta_branched.gzipped
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
                ch_bismark_index_branched.gzipped
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
                ch_bwameth_index_branched.gzipped
            )

            ch_bwameth_index = ch_bwameth_index_branched.unzipped.mix(UNTAR.out.untar)
            ch_versions      = ch_versions.mix(UNTAR.out.versions)
        } else {
            if (use_mem2) {
                BWAMETH_INDEX (
                    ch_fasta,
                    true
                )
            } else {
                BWAMETH_INDEX (
                    ch_fasta,
                    false
                )
            }
            ch_bwameth_index = BWAMETH_INDEX.out.index
            ch_versions      = ch_versions.mix(BWAMETH_INDEX.out.versions)
        }
    }

    // TODO: Here we could check if ch_or_val_bwamem_index is empty or not
    // if it is empty, we can run the BWA_INDEX subworkflow
    // if it is not empty, we need to validate the index (file or link)

    else if (params.aligner == 'bwamem'){
        log.info "BWA index not provided. Generating BWA index from FASTA file."
        /*
         * Generate BWA index from FASTA file
         */
        if (bwamem_index) {
            // TODO: Validate the BWA index
            ch_bwamem_index = bwamem_index //.map { meta, index ->
            //     if (index.toString().endsWith('.bwt')) {
            //         [meta, index]
            //     } else {
            //         error "BWA index file ${index} is not valid. It should end with .bwt"
            //     }
            // }
        } else {
            BWA_INDEX(
                ch_fasta
            )
            ch_bwamem_index = BWA_INDEX.out.index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    /*
    * Generate fasta index if not supplied for bwameth workflow or picard collecthsmetrics tool
    */
    if (aligner == 'bwameth' || aligner == 'bwamem' || collecthsmetrics) {
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
    bwamem_index  = ch_bwamem_index // channel: [ val(meta), [ bwamem index ] ]
    versions      = ch_versions      // channel: [ versions.yml ]
}
