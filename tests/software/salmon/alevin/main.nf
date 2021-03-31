#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_ALEVIN } from '../../../../software/salmon/alevin/main.nf' addParams( options: [:] )
include { SALMON_INDEX }  from '../../../../software/salmon/index/main.nf'  addParams( options: [args: "--gencode"] )
include { GFFREAD }       from '../../../../software/gffread/main.nf'       addParams( options: [args: "--table transcript_id,gene_id"] )


workflow test_salmon_alevin {

    input = [ [ id:'test', single_end:false ], // meta map
              [file("${launchDir}/tests/data/delete_me/salmon_alevin/testdata/S10_L001_R1_001.fastq.gz", checkIfExists: true),
               file("${launchDir}/tests/data/delete_me/salmon_alevin/testdata/S10_L001_R2_001.fastq.gz", checkIfExists: true)]
              ]


    transcriptome_fasta     = file("${launchDir}/tests/data/delete_me/salmon_alevin/gencode.v37.transcripts.fa", checkIfExists: true)
    genome_fasta            = file("${launchDir}/tests/data/delete_me/salmon_alevin/GRCh38.p13.genome.fa", checkIfExists: true)
    gtf                     = file("${launchDir}/tests/data/delete_me/salmon_alevin/gencode.v37.annotation.gtf", checkIfExists: true)

    GFFREAD( gtf )
    SALMON_INDEX ( transcriptome_fasta )
    SALMON_ALEVIN ( input, SALMON_INDEX.out.index, GFFREAD.out.gtf )
}
