// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_STATS     } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH           } from '../../../modules/nf-core/mosdepth/main'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {

    take:
    cram                          // channel: [mandatory] [ meta, cram, crai ]
    fasta                         // channel: [mandatory] [ fasta ]
    intervals

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow
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

