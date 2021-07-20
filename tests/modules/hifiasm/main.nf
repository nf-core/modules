#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HIFIASM } from '../../../modules/hifiasm/main.nf' addParams( options: [args:'-f0'] )

/* 
 * Test with long reads only
 */
workflow test_hifiasm_hifi_only {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['pacbio']['test_hifi_fastq_gz'], 
              checkIfExists: true) ] 
            ]

    HIFIASM ( input, [], [], false )
}

/* 
 * Test with parental reads for phasing
 */
workflow test_hifiasm_with_parental_reads {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['pacbio']['test_hifi_fastq_gz'], 
              checkIfExists: true) ] 
            ]
    paternal_kmer_dump = file(params.test_data['homo_sapiens']['illumina']['test_yak'])
    maternal_kmer_dump = file(params.test_data['homo_sapiens']['illumina']['test2_yak'])
    HIFIASM ( input, paternal_kmer_dump, maternal_kmer_dump, true )
}
