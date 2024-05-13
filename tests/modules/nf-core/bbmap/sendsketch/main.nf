#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_SENDSKETCH } from '../../../../../modules/nf-core/bbmap/sendsketch/main.nf'

workflow test_bbmap_sendsketch_single {
    
    input = [  [ id:'test', single_end:true], // meta map 
               [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
               ]
            ]

    BBMAP_SENDSKETCH ( input )
}

workflow test_bbmap_sendsketch_paired {

    input = [  [ id:'test', single_end:false], // meta map 
               [ file(params.test_data['homo_sapiens']['illumina']['test2_1_fastq_gz'], checkIfExists: true),
                 file(params.test_data['homo_sapiens']['illumina']['test2_2_fastq_gz'], checkIfExists: true)
               ]
            ]

    BBMAP_SENDSKETCH ( input )
}

