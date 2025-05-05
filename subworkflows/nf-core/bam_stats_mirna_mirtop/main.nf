include { MIRTOP_GFF        } from '../../../modules/nf-core/mirtop/gff'
include { MIRTOP_COUNTS     } from '../../../modules/nf-core/mirtop/counts'
include { MIRTOP_EXPORT     } from '../../../modules/nf-core/mirtop/export'
include { MIRTOP_STATS      } from '../../../modules/nf-core/mirtop/stats'


workflow BAM_STATS_MIRNA_MIRTOP {

    take:
    ch_bam          // channel: [ val(meta), [ bam ] ]
    ch_hairpin      // channel: [ val(meta), [ hairpin ] ]
    ch_gtf_species  // channel: [ val(meta), [ gtf ], val(species) ]

    main:

    ch_versions = Channel.empty()

    MIRTOP_GFF ( ch_bam, ch_hairpin, ch_gtf_species )
    ch_versions = ch_versions.mix(MIRTOP_GFF.out.versions)

    MIRTOP_COUNTS ( MIRTOP_GFF.out.gff, ch_hairpin, ch_gtf_species )
    ch_versions = ch_versions.mix(MIRTOP_COUNTS.out.versions)

    MIRTOP_EXPORT ( MIRTOP_GFF.out.gff, ch_hairpin, ch_gtf_species )
    ch_versions = ch_versions.mix(MIRTOP_EXPORT.out.versions)

    MIRTOP_STATS ( MIRTOP_GFF.out.gff )
    ch_versions = ch_versions.mix(MIRTOP_STATS.out.versions)

    emit:
    isomirs        = MIRTOP_EXPORT.out.tsv                  // channel: [ val(meta), [ tsv ] ]
    counts         = MIRTOP_COUNTS.out.tsv                  // channel: [ val(meta), [ tsv ] ]
    stats_txt      = MIRTOP_STATS.out.txt                   // channel: [ val(meta), [ txt ] ]
    stats_log      = MIRTOP_STATS.out.log                   // channel: [ val(meta), [ log ] ]
    versions       = ch_versions                            // channel: [ versions.yml ]
}

