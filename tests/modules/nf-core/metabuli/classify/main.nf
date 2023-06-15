#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SEQ as SEQTK_SEQ_PE } from '../../../../../modules/nf-core/seqtk/seq/main.nf'
include { SEQTK_SEQ as SEQTK_SEQ_PE_RV } from '../../../../../modules/nf-core/seqtk/seq/main.nf'
include { SEQTK_SEQ as SEQTK_SEQ_SE } from '../../../../../modules/nf-core/seqtk/seq/main.nf'
include { UNTAR } from '../../../../../modules/nf-core/untar/main.nf'
include { METABULI_CLASSIFY as METABULI_CLASSIFY_PE } from '../../../../../modules/nf-core/metabuli/classify/main.nf'
include { METABULI_CLASSIFY as METABULI_CLASSIFY_SE } from '../../../../../modules/nf-core/metabuli/classify/main.nf'
include { BUILD_ACC2TAXID, CREATE_TAXONOMY_FOLDER } from '../build'
include { METABULI_BUILD } from '../../../../../modules/nf-core/metabuli/build/main.nf'

workflow create_db {
  genome = file("${params.test_data_base}/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
  dmp_files = [
      file("${params.test_data_base}/delete_me/metabuli/names.dmp"),
      file("${params.test_data_base}/delete_me/metabuli/nodes.dmp")
  ]
  acc2taxid = BUILD_ACC2TAXID(genome)
  tax = CREATE_TAXONOMY_FOLDER(dmp_files) 
  METABULI_BUILD ( genome, acc2taxid, tax )

  emit:
    METABULI_BUILD.out
}



// test with single end data
workflow test_metabuli_classify_se {
    
    input = [
        [ id:'test_se', single_end:true ], // meta map
        [
          file("${params.test_data_base}/data/genomics/sarscov2/nanopore/fastq/test_2.fastq.gz", checkIfExists: true),
        ]
    ]
    
    db_archive =  create_db().out.db
    UNTAR( [[:], db_archive])
    SEQTK_SEQ_SE(input)
    METABULI_CLASSIFY_SE ( SEQTK_SEQ_SE.out.fastx , UNTAR.out.untar.map{it[1]})
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
    

    db_archive =  create_db().out.db

    UNTAR([[:], db_archive])
    //transform pe reads to fasta prior to classification
    
    input.map{meta, reads -> [meta, reads[0]]}
      .set{fw_reads}

    input.map{meta, reads -> [meta, reads[1]]}
      .set{rv_reads}

    SEQTK_SEQ_PE(fw_reads)

    SEQTK_SEQ_PE_RV(rv_reads)
    
    fastas = SEQTK_SEQ_PE.out.fastx
      .combine(SEQTK_SEQ_PE_RV.out.fastx, by: 0)
      .map{meta, read1, read2 -> [meta, [read1, read2]]}

    METABULI_CLASSIFY_PE ( fastas , UNTAR.out.untar.map{it[1]})
}
