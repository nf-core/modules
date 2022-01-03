#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_INDEX } from '../../../../modules/bbmap/bbsplit/main.nf'
include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_SPLIT } from '../../../../modules/bbmap/bbsplit/main.nf'

workflow test_bbmap_bbsplit {

    input = [ 
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    bbsplit_fasta_list = [ 
        ['human'],
        file('https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/chr22_23800000-23980000.fa', checkIfExists: true) 
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
	
    BBMAP_BBSPLIT_INDEX ( [ [:], [] ], [], fasta, bbsplit_fasta_list, true )
    BBMAP_BBSPLIT_SPLIT ( input, BBMAP_BBSPLIT_INDEX.out.index, fasta, bbsplit_fasta_list, true )
}
