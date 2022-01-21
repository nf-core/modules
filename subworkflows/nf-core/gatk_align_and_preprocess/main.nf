//
// Performs GATK best practice alignment and pre-processing of reads using BWA, GATK mergebamalignments (where necessary), markduplicates, sortsam, samtools index and BQSR.
// BWA index created from fasta file if not already provided
//
include { BWAMEM2_INDEX } from '../../../modules/bwamem2/index/main.nf'
include { SAMTOOLS_FASTQ } from '../../../modules/samtools/fastq/main.nf'
include { BWAMEM2_MEM   } from '../../../modules/bwamem2/mem/main.nf'
include { GATK4_MERGEBAMALIGNMENT } from '../../../modules/gatk4/mergebamalignment/main.nf'
include { PICARD_MARKDUPLICATES } from '../../../modules/picard/markduplicates/main.nf'
include { PICARD_SORTSAM } from '../../../modules/picard/sortsam/main.nf'
include { SAMTOOLS_INDEX } from '../../../modules/samtools/index/main.nf'
include { GATK4_BASERECALIBRATOR } from '../../../modules/gatk4/baserecalibrator/main.nf'
include { GATK4_APPLYBQSR } from '../../../modules/gatk4/applybqsr/main.nf'

workflow GATK_ALIGN_AND_PREPROCESS {
    take:
    input       // channel: [ val(meta), [ input ], intervals ]
    fasta               // channel: /path/to/reference/fasta
    fai                 // channel: /path/to/reference/fasta/index
    dict                // channel: /path/to/reference/fasta/dictionary
    bwaindex            // channel: name for panel of normals
    is_ubam
    sort_order
    knownsites
    knownsites_tbi

    main:
    ch_versions      = Channel.empty()

    //
    //If no bwa index has been provided, use bwamem2 index to create one from the fasta file
    //
    if (bwaindex) {
        ch_bwa_index = bwaindex
    } else {
        BWAMEM2_INDEX( fasta )
        ch_bwa_index = BWAMEM2_INDEX.out.index
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    }

    //
    //If the user inputs a unaligned bam file, it must be converted to fastq format, then aligned with bwamem2. Mergebamalignment is then used to add useful info from unaligned to the aligned bam.
    //If the user inputs a fastq file(s), the they are aligned with bwamem2 and mergebamalignment is skipped, since there is no unaligned bam file to extract info from.
    //
    if (is_ubam) {
        ch_ubam = channel.from(input)
        //
        //convert unaligned bam to fastq format
        //
        SAMTOOLS_FASTQ( ch_ubam )
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
        ch_bwa_in = SAMTOOLS_FASTQ.out.fastq.collect()

        //
        //Align reads using bwamem2 mem
        //
        BWAMEM2_MEM ( ch_bwa_in, ch_bwa_index, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        ch_markdup_in = BWAMEM2_MEM.out.bam.collect()
    } else {
        //
        //Align reads using bwamem2 mem
        //
        ch_bwa_in = channel.from(input)
        BWAMEM2_MEM ( ch_bwa_in, ch_bwa_index, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        ch_markdup_in = BWAMEM2_MEM.out.bam.collect()
    }


    //
    //use picard markduplicates to mark duplicates in the alignment bams
    //
    PICARD_MARKDUPLICATES ( ch_markdup_in )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    ch_sortsam_in = PICARD_MARKDUPLICATES.out.bam.collect()

    //
    //Bam files sorted using picard sortsam.
    //
    PICARD_SORTSAM ( ch_sortsam_in, sort_order )
    ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions)
    ch_samindex_in = PICARD_SORTSAM.out.bam.collect()

    //
    //Index for sorted bam file made using samtools index
    //
    SAMTOOLS_INDEX (ch_samindex_in)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_bai = SAMTOOLS_INDEX.out.bai.collect()

    //
    //Perform first pass of BQSR using gatk baserecalibrator.
    //
    ch_baserecal_in = ch_samindex_in.combine(ch_bai, by: 0)
    GATK4_BASERECALIBRATOR( ch_baserecal_in, fasta, fai, dict, knownsites, knownsites_tbi )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)
    ch_bqsrtable = GATK4_BASERECALIBRATOR.out.table.collect()

    //
    //Perform second pass of BQSR using gatk applybqsr.
    //
    ch_bqsr_in = ch_baserecal_in.combine(ch_bqsrtable, by: 0)
    GATK4_APPLYBQSR( ch_bqsr_in, fasta, fai, dict )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    emit:
    bwa_index_out       = ch_bwa_index

    bwa_mem_out         = BWAMEM2_MEM.out.bam.collect()                     // channel: [ val(meta), [ vcf ] ]
    // mutect2_index    = GATK4_MUTECT2.out.tbi.collect()                     // channel: [ val(meta), [ tbi ] ]
    // mutect2_stats    = GATK4_MUTECT2.out.stats.collect()                   // channel: [ val(meta), [ stats ] ]

    // genomicsdb       = GATK4_GENOMICSDBIMPORT.out.genomicsdb               // channel: [ val(meta), [ genomicsdb ] ]

    // pon_vcf          = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf           // channel: [ val(meta), [ vcf.gz ] ]
    // pon_index        = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi           // channel: [ val(meta), [ tbi ] ]

    // versions         = ch_versions                                         // channel: [ versions.yml ]
}
