#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPADES } from '../../../software/spades/main.nf' addParams( spades_hmm: false ,options: [:] )

workflow test_spades_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/generic/fastq/test_R1.fastq.gz", checkIfExists: true) ] 
            ]
    coronaspades = false

    SPADES ( input, [], coronaspades )
}

workflow test_spades_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/generic/fastq/test_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/generic/fastq/test_R2.fastq.gz", checkIfExists: true) ] 
            ]
    coronaspades = false

    SPADES ( input, [], coronaspades )
}

workflow test_coronospades_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true) ] 
            ]
    coronaspades = true

    SPADES ( input, [], coronaspades )
}

workflow test_coronospades_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ] 
            ]
    coronaspades = true

    SPADES ( input, [], coronaspades )
}

