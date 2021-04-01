#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOSDEPTH } from '../../../software/mosdepth/main.nf' addParams( options: [:] )

workflow test_mosdepth {
    input  = [ [ id:'test', single_end:true ],
               [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) ],
               [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) ] 
             ]
    dummy  = file("dummy_file.txt")
    window = 100

    MOSDEPTH ( input, dummy, window )
}
