#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNICYCLER } from '../../../software/unicycler/main.nf' addParams( options: [:] )

workflow test_unicycler_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/generic/fastq/test_R1.fastq.gz", checkIfExists: true) ] 
            ]

    UNICYCLER ( input )
}

workflow test_unicycler_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/generic/fastq/test_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/generic/fastq/test_R2.fastq.gz", checkIfExists: true) ] 
            ]

    UNICYCLER ( input )
}
