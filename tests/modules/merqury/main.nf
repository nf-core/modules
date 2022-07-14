#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERYL_COUNT    } from '../../../modules/meryl/count/main.nf'
include { MERYL_UNIONSUM } from '../../../modules/meryl/unionsum/main.nf'
include { MERQURY        } from '../../../modules/merqury/main.nf'

workflow test_merqury {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    assembly = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MERYL_COUNT ( input )
    MERYL_UNIONSUM ( MERYL_COUNT.out.meryl_db )
    MERQURY ( MERYL_UNIONSUM.out.meryl_db.join( Channel.value( assembly ) ) )
}
