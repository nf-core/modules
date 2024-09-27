include { BWA_INDEX         } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM           } from '../../../modules/nf-core/bwa/mem/main'
include { SAMBLASTER        } from '../../../modules/nf-core/samblaster/main'

workflow FASTQ_BWA_MEM_SAMBLASTER {

    take:
    ch_fastq                // channel: [ val(meta), [ fq ] ]; meta ~ [ id: 'sample' ]
    ch_reference            // channel: [ val(meta2), fasta, index ]; fasta | index; meta2 ~ [ id: 'genome' ]
                            // Each item from ch_fastq is combined with each item from ch_reference

    main:
    ch_versions             = Channel.empty()

    ch_has_index            = ch_reference
                            | branch { meta2, fasta, index ->
                                yes: index
                                no: !index
                            }

    // MODULE: BWA_INDEX
    BWA_INDEX ( ch_has_index.no.map { meta2, fasta, index -> [ meta2, fasta ] } )

    ch_bwa_index            = BWA_INDEX.out.index
                            | mix(
                                ch_has_index.yes
                                | map { meta2, fasta, index ->
                                    [ meta2, index ]
                                }
                            )

    ch_versions             = ch_versions.mix(BWA_INDEX.out.versions.first())

    // MODULE: BWA_MEM
    ch_mem_inputs           = ch_fastq
                            | combine(
                                ch_bwa_index
                            )
                            | map { meta, fq, meta2, index ->
                                [ meta + [ ref_id: meta2.id ], fq, index ]
                            }

    def sort_bam            = false
    BWA_MEM(
        ch_mem_inputs.map { meta, fq, index -> [ meta, fq ] },
        ch_mem_inputs.map { meta, fq, index -> [ [], index ] },
        [ [], [] ],
        sort_bam
    )

    ch_mem_bam              = BWA_MEM.out.bam
    ch_versions             = ch_versions.mix(BWA_MEM.out.versions.first())

    // MODULE: SAMBLASTER
    SAMBLASTER ( ch_mem_bam )

    ch_blasted_bam          = SAMBLASTER.out.bam
    ch_versions             = ch_versions.mix(SAMBLASTER.out.versions.first())

    emit:
    bam                     = SAMBLASTER.out.bam    // channel: [ val(meta), bam ]; meta ~ [ id: 'sample', ref_id: 'genome' ]
    versions                = ch_versions           // channel: [ versions.yml ]
}
