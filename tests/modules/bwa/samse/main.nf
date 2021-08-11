#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf' addParams( options: [:] )
include { BWA_ALN } from '../../../../modules/bwa/aln/main.nf' addParams( options: [:] )
include { BWA_SAMSE } from '../../../../modules/bwa/samse/main.nf' addParams( options: [:] )

workflow test_bwa_samse {

    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX ( fasta )
    BWA_ALN ( input, BWA_INDEX.out.index )
    BWA_SAMSE ( BWA_ALN.out.sai, BWA_INDEX.out.index )
}
