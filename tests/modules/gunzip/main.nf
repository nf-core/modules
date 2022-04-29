#!/usr/bin/env nextflow



include { GUNZIP } from '../../../modules/gunzip/main.nf'

workflow test_gunzip {
    input = [ [],
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]

    GUNZIP ( input )
}
