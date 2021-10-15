#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHYLOFLASH } from '../../../modules/phyloflash/main.nf' addParams( options: ['args': ''] )

process STUB_PHYLOFLASH_DATABASE {
    output:
    path("ref")                   , emit: silva_database
    path("UniVec")                , emit: univec_database

    stub:
    """
    mkdir ref
    touch UniVec
    """
}


workflow test_phyloflash_single_end {

    STUB_PHYLOFLASH_DATABASE()
    
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]]

    PHYLOFLASH ( input, STUB_PHYLOFLASH_DATABASE.out.silva_database,  STUB_PHYLOFLASH_DATABASE.out.univec_database)
}



workflow test_phyloflash {

    STUB_PHYLOFLASH_DATABASE()

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    PHYLOFLASH ( input, STUB_PHYLOFLASH_DATABASE.out.silva_database,  STUB_PHYLOFLASH_DATABASE.out.univec_database)
}
