//
// Alignment with bwa aln and sort
//

include { BWA_ALN            } from '../../../modules/nf-core/bwa/aln/main'
include { BWA_SAMSE          } from '../../../modules/nf-core/bwa/samse/main'
include { BWA_SAMPE          } from '../../../modules/nf-core/bwa/sampe/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_BWAALN {

    take:
    ch_reads // channel (mandatory): [ val(meta), path(reads) ]. subworkImportant: meta REQUIRES single_end` entry!
    ch_index // channel (mandatory): [ val(meta), path(index) ]

    main:

    ch_versions = Channel.empty()

    BWA_ALN ( ch_reads, ch_index )
    ch_versions = ch_versions.mix( BWA_ALN.out.versions.first() )

    ch_sai_for_bam = ch_reads
                        .join ( BWA_ALN.out.sai )
                        .branch {
                            meta, reads, sai ->
                                pe: !meta.single_end
                                se: meta.single_end
                        }

    BWA_SAMPE ( ch_sai_for_bam.pe, ch_index )
    ch_versions = ch_versions.mix( BWA_SAMPE.out.versions.first() )

    BWA_SAMSE ( ch_sai_for_bam.se, ch_index )
    ch_versions = ch_versions.mix( BWA_SAMSE.out.versions.first() )

    ch_bam_for_index = BWA_SAMPE.out.bam.mix( BWA_SAMSE.out.bam )

    SAMTOOLS_INDEX ( ch_bam_for_index )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = ch_bam_for_index           // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai     // channel: [ val(meta), path(bai) ]
    csi      = SAMTOOLS_INDEX.out.csi     // channel: [ val(meta), path(csi) ]

    versions = ch_versions                // channel: [ path(versions.yml) ]
}

