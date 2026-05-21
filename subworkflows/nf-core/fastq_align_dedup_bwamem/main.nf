include { BAM_SORT_STATS_SAMTOOLS       } from '../../nf-core/bam_sort_stats_samtools/main'
include { FASTQ_ALIGN_BWA               } from '../../nf-core/fastq_align_bwa/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES         } from '../../../modules/nf-core/picard/markduplicates/main'
include { PARABRICKS_FQ2BAM             } from '../../../modules/nf-core/parabricks/fq2bam/main'
include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'
workflow FASTQ_ALIGN_DEDUP_BWAMEM {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]
    ch_fasta_fai // channel: [ val(meta), fasta, fai ]
    ch_bwamem_index // channel: [ val(meta), [ bwamem index ] ]
    skip_deduplication // boolean: whether to deduplicate alignments
    use_gpu // boolean: whether to use GPU accelerated alignment
    output_fmt // string: output format for parabricks fq2bam (e.g., 'bam' or 'cram')
    interval_file // channel: [ val(meta), [ interval file ] ]
    known_sites // channel: [ val(meta), [ known sites ] ]

    main:
    ch_alignment = channel.empty()
    ch_alignment_index = channel.empty()
    ch_flagstat = channel.empty()
    ch_stats = channel.empty()
    ch_idxstats = channel.empty()
    ch_picard_metrics = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] }
    /*
    Align with parabricks GPU enabled fq2bam implementation of bwa-mem
    */
    if (use_gpu) {
        PARABRICKS_FQ2BAM(
            ch_reads,
            ch_fasta,
            ch_bwamem_index,
            interval_file,
            known_sites,
            output_fmt,
        )
        ch_alignment = PARABRICKS_FQ2BAM.out.bam
        BAM_SORT_STATS_SAMTOOLS(
            ch_alignment,
            ch_fasta_fai,
        )
        ch_alignment = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_alignment_index = BAM_SORT_STATS_SAMTOOLS.out.index
        ch_stats = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    }
    else {
        FASTQ_ALIGN_BWA(
            ch_reads,
            ch_bwamem_index,
            true,
            ch_fasta_fai
        )
        ch_alignment = FASTQ_ALIGN_BWA.out.bam
    }
    if (!skip_deduplication) {
        /*
         * Run Picard AddOrReplaceReadGroups to add read group (RG) to reads in bam file
         */        PICARD_ADDORREPLACEREADGROUPS(
            ch_alignment,
            ch_fasta_fai,
        )
        /*
         * Run Picard MarkDuplicates to mark duplicates
         */        PICARD_MARKDUPLICATES(
            PICARD_ADDORREPLACEREADGROUPS.out.bam,
            ch_fasta_fai,
        )
        /*
         * Run samtools index on deduplicated alignment
         */        SAMTOOLS_INDEX(
            PICARD_MARKDUPLICATES.out.bam
        )
        ch_alignment = PICARD_MARKDUPLICATES.out.bam
        ch_alignment_index = SAMTOOLS_INDEX.out.index
        ch_picard_metrics = PICARD_MARKDUPLICATES.out.metrics
    }
    /*
     * Collect MultiQC inputs
     */    ch_multiqc_files = ch_picard_metrics
        .collect { _meta, metrics -> metrics }
        .mix(ch_flagstat.collect { _meta, flagstat -> flagstat })
        .mix(ch_stats.collect { _meta, stats -> stats })
        .mix(ch_idxstats.collect { _meta, stats -> stats })

    emit:
    bam               = ch_alignment // channel: [ val(meta), [ bam ]       ]
    index             = ch_alignment_index // channel: [ val(meta), [ index ]     ]
    samtools_flagstat = ch_flagstat // channel: [ val(meta), [ flagstat ]  ]
    samtools_stats    = ch_stats // channel: [ val(meta), [ stats ]     ]
    samtools_idxstats = ch_idxstats // channel: [ val(meta), [ idxstats ]  ]
    picard_metrics    = ch_picard_metrics // channel: [ val(meta), [ metrics ]   ]
    multiqc           = ch_multiqc_files // channel: [ *{html,txt}              ]
}
