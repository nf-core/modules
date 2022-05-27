// runs SAMPLE_QC from either reads or bam files
// or both after alignment

include { SAMTOOLS_SORT     } from '../../../modules/samtools/sort/main'
include { SAMTOOLS_INDEX    } from '../../../modules/samtools/index/main'
include { SAMTOOLS_STATS    } from '../../../modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/samtools/flagstat/main'
include { QUALIMAP_BAMQC    } from '../../../modules/qualimap/bamqc/main'



workflow BAM_QC {

    take:
    bam        // channel: [mandatory] [ val(meta), path(bam) ]
    fasta      // channel: [mandatory] path(fasta)
    gff        // channel: [optional] path(gff)

    main:
    ch_versions = Channel.empty()

    // samtools stats block needs the bam file to be sorted
    // and indexed

    SAMTOOLS_SORT ( bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // additionally, the modules require a single channel containing
    // both the bam file and its index
    // which is the reason for creating an additional channel that joins
    // both of them:

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .set { bam_bai }

    SAMTOOLS_STATS ( bam_bai, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    SAMTOOLS_FLAGSTAT ( bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_IDXSTATS ( bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    // qualimap requires the original bam file
    // but also a GFF file with the regions to run the QC on
    QUALIMAP_BAMQC ( bam, gff )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    qualimap = QUALIMAP_BAMQC.out.results     // channel: [ val(meta), [ results ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
