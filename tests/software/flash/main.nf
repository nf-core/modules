#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FLASH } from '../../../software/flash/main.nf' addParams( options: [args:'-m 20 -M 100'] )

workflow test_flash {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ] 
            ]

    FLASH ( input )
}
