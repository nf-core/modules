#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { YARA_INDEX } from '../../../../software/yara/index/main.nf' addParams( options: ['args': '-e 3'] )
include { YARA_MAPPER } from '../../../../software/yara/mapper/main.nf' addParams( options: ['args': '-e 3'] )

workflow test_yara_single_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    YARA_INDEX ( fasta )

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true) ]

    YARA_MAPPER ( input, YARA_INDEX.out.index )
}

workflow test_yara_paired_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    YARA_INDEX ( fasta )

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ] ]

    YARA_MAPPER ( input, YARA_INDEX.out.index )
}
