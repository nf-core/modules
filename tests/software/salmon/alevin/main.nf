#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_ALEVIN } from '../../../../software/salmon/alevin/main.nf' addParams( options: [:] )
include { GFFREAD }       from '../../../../software/gffread/main.nf'       addParams( options: [args: "--table transcript_id,gene_id"] )


workflow test_salmon_alevin {

    input = [ [ id:'test', single_end:false ], // meta map
              [file("${launchDir}/tests/data/delete_me/salmon_alevin_data/S10_L001_R1_001.fastq.gz", checkIfExists: true),
               file("${launchDir}/tests/data/delete_me/salmon_alevin_data/S10_L001_R2_001.fastq.gz", checkIfExists: true)]
              ]


    transcriptome_fasta     = file("${launchDir}/tests/data/delete_me/salmon_alevin_data/gencode.vM26.transcripts.fa", checkIfExists: true)
    gtf                     = file("${launchDir}/tests/data/delete_me/salmon_alevin_data/gencode.vM26.annotation.gtf", checkIfExists: true)
    salmon_index            = file("${launchDir}/tests/data/delete_me/salmon_alevin_data/salmon_index", checkIfExists: true)
    protocol                = "chromium"


    GFFREAD( gtf )
    SALMON_ALEVIN ( input, salmon_index, GFFREAD.out.gtf, protocol )
}
