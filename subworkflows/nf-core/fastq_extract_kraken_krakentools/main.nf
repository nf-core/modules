include { KRAKEN2_KRAKEN2                } from '../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS } from '../../../modules/nf-core/krakentools/extractkrakenreads/main'

workflow FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS {

    take:
    ch_reads  // channel: [ val(meta), path(reads) ]
    ch_db     // channel: [ path db ]
    val_taxid // string: taxonomic ids, separated by spaces

    main:
    ch_versions = Channel.empty()

    KRAKEN2_KRAKEN2 ( ch_reads, ch_db, true, true )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    KRAKENTOOLS_EXTRACTKRAKENREADS ( val_taxid, KRAKEN2_KRAKEN2.out.classified_reads_assignment, KRAKEN2_KRAKEN2.out.classified_reads_fastq, KRAKEN2_KRAKEN2.out.report )
    ch_versions = ch_versions.mix( KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions.first() )

    emit:
    kraken2_report          = KRAKEN2_KRAKEN2.out.report                                 // channel: [ val(meta), path ]
    extracted_kraken2_reads = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads // channel: [ val(meta), [ fastq.gz/fasta.gz ] ]
    multiqc_files           = KRAKEN2_KRAKEN2.out.report.map{it[1]}                      // channel: [ path ]
    versions                = ch_versions                                                // channel: [ versions.yml ]
}
