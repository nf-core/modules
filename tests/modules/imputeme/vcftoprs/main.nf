#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IMPUTEME_VCFTOPRS } from '../../../../modules/imputeme/vcftoprs/main.nf' addParams( options: [:] )

workflow test_imputeme_vcftoprs {
    
    input = [ [ id:'test'], // meta map
              file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true) ]

    IMPUTEME_VCFTOPRS ( input )
}
