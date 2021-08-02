#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf' addParams( options: [:] )
include { BWA_ALN } from '../../../../modules/bwa/aln/main.nf' addParams( options: [:] )
include { BWA_SAMPE } from '../../../../modules/bwa/sampe/main.nf' addParams( options: [:] )

workflow test_bwa_sampe {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX ( fasta )
    BWA_ALN ( input, BWA_INDEX.out.index )
    BWA_SAMPE ( BWA_ALN.out.sai, BWA_INDEX.out.index )
}
