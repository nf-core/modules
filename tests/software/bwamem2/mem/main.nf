#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../../../software/bwamem2/index/main.nf' addParams( options: [:] )
include { BWAMEM2_MEM   } from '../../../../software/bwamem2/mem/main.nf'   addParams( options: [:] )

/*
 * Test with single-end data
 */
workflow test_bwamem2_mem_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true) ] 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BWAMEM2_INDEX ( fasta )
    BWAMEM2_MEM ( input, BWAMEM2_INDEX.out.index )
}

/*
 * Test with paired-end data
 */
workflow test_bwamem2_mem_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ] 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BWAMEM2_INDEX ( fasta )
    BWAMEM2_MEM ( input, BWAMEM2_INDEX.out.index )
}
