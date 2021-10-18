#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_LIFTOVER } from '../../../../modules/ucsc/liftover/main.nf' addParams( options: [:] )

workflow test_ucsc_liftover {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)]
    chain =  file("http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz")

    UCSC_LIFTOVER ( input, chain )
}
