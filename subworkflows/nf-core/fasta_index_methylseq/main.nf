include { UNTAR as UNTAR_BISMARK    } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BWAMETH    } from '../../../modules/nf-core/untar/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION as BISMARK_GENOMEPREPARATION_BOWTIE } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BISMARK_GENOMEPREPARATION as BISMARK_GENOMEPREPARATION_HISAT } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX             } from '../../../modules/nf-core/bwameth/index/main'
include { BWA_INDEX                 } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'

workflow FASTA_INDEX_METHYLSEQ {

    take:
    fasta            // channel: [ val(meta), [ fasta ] ]
    fasta_index      // channel: [ val(meta), [ fasta index ] ]
    bismark_index    // channel: [ val(meta), [ bismark index ] ]
    bwameth_index    // channel: [ val(meta), [ bwameth index ] ]
    bwamem_index     // channel: [ val(meta), [ bwamem index ] ]
    aligner          // string: bismark, bismark_hisat, bwameth or bwamem
    collecthsmetrics // boolean: whether to run picard collecthsmetrics
    use_mem2         // boolean: generate mem2 index if no index provided, and bwameth is selected

    main:

    ch_fasta         = channel.empty()
    ch_fasta_index   = channel.empty()
    ch_bismark_index = channel.empty()
    ch_bwameth_index = channel.empty()
    ch_bwamem_index  = channel.empty()
    ch_versions      = channel.empty()

    // Check if fasta file is gzipped and decompress if needed
    fasta
        .branch { _meta, file ->
            gzipped: file.toString().endsWith('.gz')
            unzipped: true
        }
        .set { ch_fasta_branched }

    GUNZIP (
        ch_fasta_branched.gzipped
    )

    ch_fasta    = ch_fasta_branched.unzipped.mix(GUNZIP.out.gunzip)

    // Aligner: bismark or bismark_hisat
    if( aligner =~ /bismark/ ){
        /*
         * Generate bismark index if not supplied
         */
        if (bismark_index) {
            // Handle channel-based bismark index
            bismark_index
                .branch { _meta, file ->
                    gzipped: file.toString().endsWith('.gz')
                    unzipped: true
                }
                .set { ch_bismark_index_branched }

            UNTAR_BISMARK (
                ch_bismark_index_branched.gzipped
            )

            ch_bismark_index = ch_bismark_index_branched.unzipped.mix(UNTAR_BISMARK.out.untar)
        } else {

            if( aligner == "bismark_hisat") {
                BISMARK_GENOMEPREPARATION_HISAT (
                    ch_fasta
                )
                ch_bismark_index = BISMARK_GENOMEPREPARATION_HISAT.out.index
                ch_versions      = ch_versions.mix(BISMARK_GENOMEPREPARATION_HISAT.out.versions)
            } else {
                BISMARK_GENOMEPREPARATION_BOWTIE (
                    ch_fasta
                )
                ch_bismark_index = BISMARK_GENOMEPREPARATION_BOWTIE.out.index
                ch_versions      = ch_versions.mix(BISMARK_GENOMEPREPARATION_BOWTIE.out.versions)
            }
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
                .branch { _meta, file ->
                    gzipped: file.toString().endsWith('.gz')
                    unzipped: true
                }
                .set { ch_bwameth_index_branched }

            UNTAR_BWAMETH (
                ch_bwameth_index_branched.gzipped
            )

            ch_bwameth_index = ch_bwameth_index_branched.unzipped.mix(UNTAR_BWAMETH.out.untar)
        } else {
            BWAMETH_INDEX (
                ch_fasta,
                use_mem2
            )
            ch_bwameth_index = BWAMETH_INDEX.out.index
            ch_versions      = ch_versions.mix(BWAMETH_INDEX.out.versions)
        }
    }


    else if ( aligner == 'bwamem' ){
        /*
         * Generate BWA index from FASTA file
         */
        if (bwamem_index) {
            // Handle channel-based bwamem index
            bwamem_index
                .branch { _meta, file ->
                    gzipped: file.toString().endsWith('.gz')
                    unzipped: true
                }
                .set { ch_bwamem_index_branched }

            UNTAR_BISMARK (
                ch_bwamem_index_branched.gzipped
            )

            ch_bwamem_index = ch_bwamem_index_branched.unzipped.mix(UNTAR_BISMARK.out.untar)
        } else {
            log.info "BWA index not provided. Generating BWA index from FASTA file."
            BWA_INDEX (
                ch_fasta
            )
            ch_bwamem_index = BWA_INDEX.out.index
        }
    }

    /*
    * Generate fasta index if not supplied for bwameth workflow or picard collecthsmetrics tool
    */
    if (aligner == 'bwameth' || aligner == 'bwamem' || collecthsmetrics) {
        // already existing fasta index
        if (fasta_index) {
            ch_fasta_index = fasta_index
        } else {
            log.info "Fasta index not provided. Generating fasta index from FASTA file."
            SAMTOOLS_FAIDX (
                ch_fasta.combine(channel.of([[]])),
                false
            )
            ch_fasta_index = SAMTOOLS_FAIDX.out.fai
            // samtools/faidx version emitted into the topic channel
        }
    }

    emit:
    fasta         = ch_fasta         // channel: [ val(meta), [ fasta ] ]
    fasta_index   = ch_fasta_index   // channel: [ val(meta), [ fasta index ] ]
    bismark_index = ch_bismark_index // channel: [ val(meta), [ bismark index ] ]
    bwameth_index = ch_bwameth_index // channel: [ val(meta), [ bwameth index ] ]
    bwamem_index  = ch_bwamem_index  // channel: [ val(meta), [ bwamem index ] ]
    versions      = ch_versions      // channel: [ versions.yml ]
}
