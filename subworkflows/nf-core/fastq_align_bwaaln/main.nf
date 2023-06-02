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

    // Ensure multiple references and read combinations are correctly associated throughout the subworkflow
    ch_prepped_input = ch_reads
                        .combine(ch_index)
                        .map{
                            meta, reads, meta_index, index ->

                            def key_read_ref = [ meta.id.join(meta_index.id, "_")]

                            [ meta + [key_read_ref: key_read_ref], reads, meta_index + [key_read_ref: key_read_ref], index  ]
                        }
                        .branch {
                            meta, reads, meta_index, index ->
                            reads: [ meta      , reads ]
                            index: [ meta_index, index ]
                        }

    // Alignment and conversion to bAM
    BWA_ALN ( ch_prepped_input.reads, ch_prepped_input.index )
    ch_versions = ch_versions.mix( BWA_ALN.out.versions.first() )

    ch_sai_for_bam = ch_prepped_input
                        .join ( BWA_ALN.out.sai )
                        .branch {
                            meta, reads, sai ->
                                pe: !meta.single_end
                                se: meta.single_end
                        }

    // ch_sai_for_bam_pe =  ch_sai_for_bam.pe
    //                         .combine{
    //                             ch_prepped_input.index
    //                         }
    //                         .dump(tag: "hello")
    //                         .mulimap {
    //                             meta, reads, meta_index, index ->
    //                                 reads: [ meta      , reads ]
    //                                 index: [ meta_index, index ]
    //                         }

    // ch_sai_for_bam_se =  ch_sai_for_bam.se
    //                         .combine{
    //                             ch_prepped_input.index
    //                         }
    //                         .mulimap {
    //                             meta, reads, meta_index, index ->
    //                                 reads: [ meta      , reads ]
    //                                 index: [ meta_index, index ]
    //                         }

    // BWA_SAMPE ( ch_sai_for_bam_pe.reads, ch_sai_for_bam_pe.index )
    // ch_versions = ch_versions.mix( BWA_SAMPE.out.versions.first() )

    // BWA_SAMSE ( ch_sai_for_bam_se.reads, ch_sai_for_bam_se.index )
    // ch_versions = ch_versions.mix( BWA_SAMSE.out.versions.first() )

    // ch_bam_for_index = BWA_SAMPE.out.bam.mix( BWA_SAMSE.out.bam )

    // SAMTOOLS_INDEX ( ch_bam_for_index )
    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // bam      = ch_bam_for_index.dump(tag: "bam")           // channel: [ val(meta), path(bam) ]
    // bai      = SAMTOOLS_INDEX.out.bai.dump(tag: "bai")     // channel: [ val(meta), path(bai) ]
    // csi      = SAMTOOLS_INDEX.out.csi.dump(tag: "csi")     // channel: [ val(meta), path(csi) ]

    versions = ch_versions                // channel: [ path(versions.yml) ]
}

