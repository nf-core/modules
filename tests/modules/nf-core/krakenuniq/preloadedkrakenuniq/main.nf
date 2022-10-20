#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKENUNIQ_DOWNLOAD } from '../../../../../modules/nf-core/krakenuniq/download/main.nf'
include { KRAKENUNIQ_BUILD    } from '../../../../../modules/nf-core/krakenuniq/build/main.nf'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ } from '../../../../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main.nf'

workflow test_krakenuniq_preloadedkrakenuniq {
    ch_fastas      = Channel.fromPath(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ch_seqid2mapid = Channel.fromPath(params.test_data['sarscov2']['metagenome']['seqid2taxid_map'], checkIfExists: true)
    KRAKENUNIQ_DOWNLOAD ( 'taxonomy' )

    ch_input = ch_fastas
                    .combine(KRAKENUNIQ_DOWNLOAD.out.output)
                    .combine(ch_seqid2mapid)
                    .map {
                        fna, tax, map ->

                        [ [id: "customdb"], fna, tax, map ]
                    }

    KRAKENUNIQ_BUILD ( ch_input )

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz', 'test_2_fastq_gz'], checkIfExists: true)
    ]

    KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( input,KRAKENUNIQ_BUILD.out.db.map{it[1]},true,true,true )
}
