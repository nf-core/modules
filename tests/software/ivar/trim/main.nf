#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_TRIM } from '../../../../software/ivar/trim/main.nf' addParams([:])

workflow test_ivar_trim {
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) 
            ]
    bed   = file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)

    IVAR_TRIM ( input, bed_file )
}
