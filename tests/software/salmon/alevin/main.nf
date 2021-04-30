#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_ALEVIN } from '../../../../software/salmon/alevin/main.nf' addParams( options: [:] )
include { GFFREAD }       from '../../../../software/gffread/main.nf'       addParams( options: [args: "--table transcript_id,gene_id"] )


workflow test_salmon_alevin {

    input = [ [ id:'test', single_end:false ], // meta map
              [file(params.test_data['homo_sapiens']['illumina']['test_1.fastq.gz'], checkIfExists: true),
               file(params.test_data['homo_sapiens']['illumina']['test_1.fastq.gz'], checkIfExists: true)]
              ]


    transcriptome_fasta     = file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)
    gtf                     = file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)
    salmon_index            = file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)
    protocol                = "chromium"
    genome_fasta            = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)


    SALMON_INDEX ( genome_fasta, transctranscriptome_fastaript_fasta )

    GFFREAD( gtf )
    SALMON_ALEVIN ( input, salmon_index, GFFREAD.out.gtf, protocol )
}
