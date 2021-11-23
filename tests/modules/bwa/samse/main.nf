#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf'
include { BWA_ALN   } from '../../../../modules/bwa/aln/main.nf'
include { BWA_SAMSE } from '../../../../modules/bwa/samse/main.nf'

workflow test_bwa_samse {

    Channel
        .fromList(
            [ [ id:'test', single_end:true ],
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
        ).collect()
        .set { input }
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX ( fasta )
    BWA_ALN ( input, BWA_INDEX.out.index )
    BWA_SAMSE ( input.join(BWA_ALN.out.sai, by:[0]), BWA_INDEX.out.index )
}
