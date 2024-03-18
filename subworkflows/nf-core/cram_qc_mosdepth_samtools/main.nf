//
// QC on CRAM
//

include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/'
include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth/'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {

    take:
    cram                          // channel: [mandatory] [ meta, cram, crai ]
    fasta                         // channel: [mandatory] [ fasta ]
    intervals

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_STATS(cram, fasta.map{ it -> [ [ id:'fasta' ], it ] })
    
    MOSDEPTH(cram.combine(intervals.map{ meta, bed -> [ bed?:[] ] }), fasta.map{ it -> [ [ id:'fasta' ], it ] })


    // Gather all reports generated
    reports = reports.mix(SAMTOOLS_STATS.out.stats)
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.regions_txt)

    // Gather versions of all tools used
    versions = versions.mix(MOSDEPTH.out.versions)
    versions = versions.mix(SAMTOOLS_STATS.out.versions.first())


    emit:
    reports

    versions = ch_versions                     // channel: [ versions.yml ]
}

