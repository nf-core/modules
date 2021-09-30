#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_INDEX } from '../../../../modules/bbmap/bbsplit/main.nf' addParams( options: [ args:'' ] )
include { BBMAP_BBSPLIT } from '../../../../modules/bbmap/bbsplit/main.nf' addParams( options: [ args:'' ] )

workflow test_bbmap_bbsplit {

    input                 = [ [ id:'test', single_end:true ], // meta map
                              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
                            ]
    fasta                 = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ch_bbsplit_fasta_list = [ ['human'],
                              file('https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/chr22_23800000-23980000.fa', checkIfExists: true) ]

    BBMAP_BBSPLIT_INDEX ( [ [:], [] ], [], fasta, ch_bbsplit_fasta_list, true )
    BBMAP_BBSPLIT ( input, BBMAP_BBSPLIT_INDEX.out.index, fasta, ch_bbsplit_fasta_list, true )
}
