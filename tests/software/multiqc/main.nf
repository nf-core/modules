#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC  } from '../../../software/fastqc/main.nf'  addParams( options: [:] )
include { MULTIQC } from '../../../software/multiqc/main.nf' addParams( options: [:] )

workflow test_multiqc {
    input = [ [ id: 'test', single_end: false ],
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ] 
            ]
    
    FASTQC  ( input )
    MULTIQC ( FASTQC.out.zip.collect { it[1] } )
}
