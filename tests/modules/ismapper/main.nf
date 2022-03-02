#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ISMAPPER } from '../../../modules/ismapper/main.nf'

workflow test_ismapper {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
        file("https://github.com/jhawkey/IS_mapper/raw/master/test/inputs/S_suis_P17.gbk", checkIfExists: true),
        file("https://github.com/jhawkey/IS_mapper/raw/master/test/inputs/ISSsu3.fasta", checkIfExists: true)
    ]

    ISMAPPER ( input )
}
