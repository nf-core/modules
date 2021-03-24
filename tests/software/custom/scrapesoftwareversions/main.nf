#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from '../../../../software/fastqc/main.nf' addParams( options: [:] )
include { CUSTOM_SCRAPESOFTWAREVERSIONS } from '../../../../software/custom/scrapesoftwareversions/main.nf' addParams( options: [:] )

workflow test_scrapesoftwareversions {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true) ]
            ]
    FASTQC ( input )
    CUSTOM_SCRAPESOFTWAREVERSIONS ( FASTQC.out.version )
}
