//
// BAM deduplication with UMI processing
//

include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME } from '../../../subworkflows/nf-core/bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME } from '../../../subworkflows/nf-core/bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../../../modules/nf-core/umitools/prepareforrsem/main'

workflow BAM_DEDUP_UMI_PROCESS {
    take:
    ch_genome_bam         // channel: [ val(meta), path(bam) ]
    ch_genome_bam_index   // channel: [ val(meta), path(bai) ]
    ch_fasta              // channel: [ path(fasta) ]
    umi_dedup_tool        // string: 'umicollapse' or 'umitools'
    umitools_dedup_stats  // boolean: whether to generate UMI-tools dedup stats
    bam_csi_index         // boolean: whether to generate CSI index
    with_transcriptome    // boolean: whether to process transcriptome BAM
    ch_transcriptome_bam  // channel: [ val(meta), path(bam) ]
    ch_transcript_fasta   // channel: [ path(transcript_fasta) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Genome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME (
            ch_genome_bam.join(ch_genome_bam_index, by: [0])
        )
        ch_dedup_genome = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME
        ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_genome.out.dedup_stats.collect{it[1]}.ifEmpty([]))
    } else if (umi_dedup_tool == "umitools") {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
            ch_genome_bam.join(ch_genome_bam_index, by: [0]),
            umitools_dedup_stats
        )
        ch_dedup_genome = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME
        ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_genome.out.deduplog.collect{it[1]})
    } else {
        error("Unknown umi_dedup_tool '${umi_dedup_tool}'")
    }

    ch_genome_bam       = ch_dedup_genome.out.bam
    ch_genome_bam_index = ch_dedup_genome.out.bai
    ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_genome.out.stats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_genome.out.flagstat.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_genome.out.idxstats.collect{it[1]})

    if (bam_csi_index) {
        ch_genome_bam_index = ch_dedup_genome.out.csi
    }
    ch_versions = ch_versions.mix(ch_dedup_genome.out.versions)

    // Initialize transcriptome-dependent channels as empty
    ch_transcriptome_bam_dedup = Channel.empty()

    // Transcriptome BAM processing
    if (with_transcriptome) {
        // Co-ordinate sort, index and run stats on transcriptome BAM
        BAM_SORT_STATS_SAMTOOLS (
            ch_transcriptome_bam,
            ch_transcript_fasta.map { [ [:], it ] }
        )
        ch_transcriptome_sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_transcriptome_sorted_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

        // Deduplicate transcriptome BAM file
        if (umi_dedup_tool == "umicollapse") {
            BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME (
                ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0])
            )
            ch_dedup_transcriptome = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME
        } else if (umi_dedup_tool == "umitools") {
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
                ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
                umitools_dedup_stats
            )
            ch_dedup_transcriptome = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME
        }

        // Name sort BAM
        SAMTOOLS_SORT (
            ch_dedup_transcriptome.out.bam,
            ch_fasta.map { [ [:], it ] }
        )

        // Branch for single-end and paired-end BAMs
        SAMTOOLS_SORT
            .out
            .bam
            .branch {
                meta, bam ->
                    single_end: meta.single_end
                        return [ meta, bam ]
                    paired_end: !meta.single_end
                        return [ meta, bam ]
            }
            .set { ch_dedup_bam }

        // Fix paired-end reads in name sorted BAM file
        UMITOOLS_PREPAREFORSALMON (
            ch_dedup_bam.paired_end.map { meta, bam -> [ meta, bam, [] ] }
        )
        ch_versions = ch_versions.mix(UMITOOLS_PREPAREFORSALMON.out.versions.first())

        ch_transcriptome_bam_dedup = ch_dedup_bam.single_end.mix(UMITOOLS_PREPAREFORSALMON.out.bam)
    }

    emit:
    genome_bam          = ch_genome_bam                  // channel: [ val(meta), path(bam) ]
    genome_bam_index    = ch_genome_bam_index            // channel: [ val(meta), path(bai) ]
    transcriptome_bam   = ch_transcriptome_bam_dedup     // channel: [ val(meta), path(bam) ]
    versions            = ch_versions                    // channel: [ path(versions.yml) ]
    multiqc_files       = ch_multiqc_files               // channel: [ path(multiqc_files) ]
}
