#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../../../modules/nf-core/untar/main.nf'
include { METABULI_CLASSIFY as METABULI_CLASSIFY_PE } from '../../../../../modules/nf-core/metabuli/classify/main.nf'
include { METABULI_CLASSIFY as METABULI_CLASSIFY_SE } from '../../../../../modules/nf-core/metabuli/classify/main.nf'

// test with single end data
workflow test_metabuli_classify_se {
    
    input = [
        [ id:'test_se', single_end:true ], // meta map
        [
          file("${params.test_data_base}/data/genomics/sarscov2/nanopore/fastq/test_2.fastq.gz", checkIfExists: true),
        ]
    ]
        
    db_archive =  file("${params.test_data_base}/data/delete_me/metabuli/metabuli_db.tar.gz", checkIfExists: true)
    UNTAR( [[:], db_archive])
    METABULI_CLASSIFY_SE ( input , UNTAR.out.untar.map{it[1]})
}

// test with paired end data
workflow test_metabuli_classify_pe {
    
    input = Channel.from(
        [[ id:'test_pe', single_end:false ], // meta map
        [
          file("${params.test_data_base}/data/genomics/sarscov2/illumina/fastq/test2_1.fastq.gz", checkIfExists: true),
          file("${params.test_data_base}/data/genomics/sarscov2/illumina/fastq/test2_2.fastq.gz", checkIfExists: true),
        ]],
        [[ id:'test_pe2', single_end:false ], // meta map
        [
          file("${params.test_data_base}/data/genomics/sarscov2/illumina/fastq/test2_1.fastq.gz", checkIfExists: true),
          file("${params.test_data_base}/data/genomics/sarscov2/illumina/fastq/test2_2.fastq.gz", checkIfExists: true),
        ]]
    )

    db_archive =  file("${params.test_data_base}/data/delete_me/metabuli/metabuli_db.tar.gz", checkIfExists: true)
    UNTAR([[:], db_archive])
    METABULI_CLASSIFY_PE ( input , UNTAR.out.untar.map{it[1]})
}
