#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHYLOFLASH } from '../../../modules/phyloflash/main.nf' addParams( options: [:] )

workflow test_phyloflash {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    PHYLOFLASH ( input )
}
