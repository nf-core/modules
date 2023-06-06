#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SEQ } from '../../../../../modules/nf-core/seqtk/seq/main.nf'
include { UNTAR } from '../../../../../modules/nf-core/untar/main.nf'
include { METABULI_CLASSIFY } from '../../../../../modules/nf-core/metabuli/classify/main.nf'

// test with single end data
workflow test_metabuli_classify_se {
    
    input = [
        [ id:'test_se', single_end:true ], // meta map
        [
          file("${params.test_data_base}/data/genomics/sarscov2/nanopore/fastq/test_2.fastq.gz", checkIfExists: true),
        ]
    ]

    db_archive = [
        [ id:'test_se'], // meta map
        file("refseq_virus.tar.gz",checkIfExists: true)
        //file("${params.test_data_base}/data/delete_me/metabuli/classify",checkIfExists: true)
    ]

    UNTAR(db_archive)
    SEQTK_SEQ(input)
    METABULI_CLASSIFY ( SEQTK_SEQ.out.fastx , UNTAR.out.untar.map{it[1]})
}

// test with paired end data
workflow test_metabuli_classify_pe {
    
    input = [
        [ id:'test_pe', single_end:false ], // meta map
        [
          file("${params.test_data_base}/data/genomics/sarscov2/illumina/fastq/test2_1.fastq.gz", checkIfExists: true),
          file("${params.test_data_base}/data/genomics/sarscov2/illumina/fastq/test2_2.fastq.gz", checkIfExists: true),
        ]
    ]

    db_archive = [
        [ id:'test_pe', single_end:false ], // meta map
        file("refseq_virus.tar.gz",checkIfExists: true)
        //file("${params.test_data_base}/data/delete_me/metabuli/classify",checkIfExists: true)
    ]

    UNTAR(db_archive)
    // TODO: transform both reads to fasta prior to classification
    //METABULI_CLASSIFY ( (fastas) , UNTAR.out.untar.map{it[1]})
}
