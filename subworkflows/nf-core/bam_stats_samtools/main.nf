nextflow.enable.types = true

//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai   // channel: records of [ meta, bam, bai ]
    ch_fasta_fai // channel: records of [ meta, fasta, fai ]

    main:
    SAMTOOLS_STATS(ch_bam_bai, ch_fasta_fai)

    SAMTOOLS_FLAGSTAT(ch_bam_bai)

    SAMTOOLS_IDXSTATS(ch_bam_bai)

    emit:
    stats    = SAMTOOLS_STATS.out    // channel: records of [ meta, stats ]
    flagstat = SAMTOOLS_FLAGSTAT.out // channel: records of [ meta, flagstat ]
    idxstats = SAMTOOLS_IDXSTATS.out // channel: records of [ meta, idxstats ]
}
