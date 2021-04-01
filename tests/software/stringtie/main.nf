#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRINGTIE as STRINGTIE_FORWARD } from '../../../software/stringtie/main.nf'  addParams( options: [ publish_dir:'test_stringtie_forward' ] )
include { STRINGTIE as STRINGTIE_REVERSE } from '../../../software/stringtie/main.nf'  addParams( options: [ publish_dir:'test_stringtie_reverse' ] )

/*
 * Test with forward strandedness
 */
workflow test_stringtie_forward {
    input = [ [ id:'test', strandedness:'forward' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) ] ]
    gtf   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true)
    
    STRINGTIE_FORWARD ( input, gtf )
}

/*
 * Test with reverse strandedness
 */
workflow test_stringtie_reverse {
    input = [ [ id:'test', strandedness:'reverse' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) ] 
            ]
    gtf   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true)
    
    STRINGTIE_REVERSE ( input, gtf )
}
