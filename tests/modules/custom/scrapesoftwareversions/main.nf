#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from '../../../../modules/fastqc/main.nf' addParams( options: [:] )
include { CUSTOM_SCRAPESOFTWAREVERSIONS } from '../../../../modules/custom/scrapesoftwareversions/main.nf' addParams( options: [:] )

workflow test_scrapesoftwareversions {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data.sarscov2.illumina.test_1_fastq_gz, checkIfExists: true) ]
            ]
    FASTQC ( input )
    CUSTOM_SCRAPESOFTWAREVERSIONS ( FASTQC.out.version )
}
