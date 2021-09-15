#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_ALIGN } from '../../../../modules/bbmap/align/main.nf' addParams( options: [:] )

workflow test_bbmap_align {
    
    input = [ [ id:'test', single_end:false ], // meta map
                [
                    file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file( params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
                ]
            ]
    ref   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_ALIGN ( input, ref )
}
