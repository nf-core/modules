//
// QC on CRAM
//

include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth'
include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    fasta_fai     // channel: [mandatory] [ meta, fasta, fai ]
    intervals
    
    main:
    reports = channel.empty()

    // Reports run on cram
    SAMTOOLS_STATS(cram, fasta_fai)

    fasta = fasta_fai.map { meta2, fasta_file, fai -> [ meta2, fasta_file ] }

    ch_mosdepth_in = cram
        .combine(intervals.map { meta, interval -> interval }.ifEmpty([ [] ]))
        .map { meta, cram_file, crai, bed -> [ meta, cram_file, crai, bed ?: [] ] }
    
    MOSDEPTH(ch_mosdepth_in, fasta, [])

    // Gather all reports generated
    reports = reports.mix(SAMTOOLS_STATS.out.stats)
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.regions_txt)

    emit:
    reports // channel: [ meta, report_file ]
}
