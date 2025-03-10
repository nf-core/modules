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

    // WARNING: You must specify in your prefix `meta.id_index` in your `modules.conf`
    // to ensure that you do not overwrite multiple BAM files from one sample mapped
    // against multiple references. This meta map is added by the subworkflow but can be removed
    // after if necessary.

    // Ensure when multiple references are provide, that reference/read combinations
    // are correctly associated throughout the subworkflow by copying the sample
    // specific metadata to the index on each combination

    ch_prepped_input = ch_reads
                        .combine(ch_index)
                        .map{
                            meta, reads, meta_index, index ->

                                // Create a combined id that includes the ids of the reads and the index used.
                                // Also keep the id of the index with a new name to avoid name collisions.
                                def key_read_ref = meta.id + "_" + meta_index.id
                                def id_index = meta_index.id

                            [ meta + [key_read_ref: key_read_ref] + [id_index: id_index], reads, meta_index + [key_read_ref: key_read_ref]  + [id_index: id_index], index  ]
                        }

    // Drop the index_meta, as the id of the index is now kept within the read meta.
    ch_preppedinput_for_bwaaln = ch_prepped_input
                        .multiMap {
                            meta, reads, meta_index, index ->
                                reads: [ meta, reads ]
                                index: [ meta, index ]
                        }


    // Set as independent channel to allow repeated joining but _with_ sample specific metadata
    // to ensure right reference goes with right sample
    ch_reads_newid = ch_prepped_input.map{ meta, reads, meta_index, index -> [ meta, reads ] }
    ch_index_newid = ch_prepped_input.map{ meta, reads, meta_index, index -> [ meta, index ] }

    // Alignment and conversion to bam
    BWA_ALN ( ch_preppedinput_for_bwaaln.reads, ch_preppedinput_for_bwaaln.index )
    ch_versions = ch_versions.mix( BWA_ALN.out.versions.first() )

    ch_sai_for_bam = ch_reads_newid
                        .join ( BWA_ALN.out.sai )
                        .branch {
                            meta, reads, sai ->
                                pe: !meta.single_end
                                se: meta.single_end
                        }

    // Split as PE/SE have different SAI -> BAM commands
    ch_sai_for_bam_pe =  ch_sai_for_bam.pe
                            .join ( ch_index_newid )
                            .multiMap {
                                meta, reads, sai, index ->
                                    reads: [ meta, reads, sai ]
                                    index: [ meta, index      ]
                            }

    ch_sai_for_bam_se =  ch_sai_for_bam.se
                            .join ( ch_index_newid )
                            .multiMap {
                                meta, reads, sai, index ->
                                    reads: [ meta, reads, sai ]
                                    index: [ meta, index      ]
                            }


    BWA_SAMPE ( ch_sai_for_bam_pe.reads, ch_sai_for_bam_pe.index )
    ch_versions = ch_versions.mix( BWA_SAMPE.out.versions.first() )

    BWA_SAMSE ( ch_sai_for_bam_se.reads, ch_sai_for_bam_se.index )
    ch_versions = ch_versions.mix( BWA_SAMSE.out.versions.first() )

    ch_bam_for_index = BWA_SAMPE.out.bam.mix( BWA_SAMSE.out.bam )

    // Index all
    SAMTOOLS_INDEX ( ch_bam_for_index )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Remove superfluous internal maps to minimise clutter as much as possible
    ch_bam_for_emit = ch_bam_for_index.map{ meta, bam -> [meta - meta.subMap('key_read_ref'), bam] }
    ch_bai_for_emit = SAMTOOLS_INDEX.out.bai.map{ meta, bai -> [meta - meta.subMap('key_read_ref'), bai] }
    ch_csi_for_emit = SAMTOOLS_INDEX.out.csi.map{ meta, csi -> [meta - meta.subMap('key_read_ref'), csi] }

    emit:
    // Note: output channels will contain meta with additional 'id_index' meta
    // value to allow association of BAM file with the meta.id of input indices
    bam      = ch_bam_for_emit     // channel: [ val(meta), path(bam) ]
    bai      = ch_bai_for_emit     // channel: [ val(meta), path(bai) ]
    csi      = ch_csi_for_emit     // channel: [ val(meta), path(csi) ]

    versions = ch_versions         // channel: [ path(versions.yml) ]
}
