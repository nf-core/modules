#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../software/gunzip/main.nf' addParams( options: [:] )

workflow test_gunzip {
    input = file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true)

    GUNZIP ( input )
}
