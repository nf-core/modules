#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_EXPORT as CNVKIT_EXPORT_BED } from '../../../../../modules/nf-core/cnvkit/export/main.nf'
include { CNVKIT_EXPORT as CNVKIT_EXPORT_VCF } from '../../../../../modules/nf-core/cnvkit/export/main.nf'
include { CNVKIT_EXPORT as CNVKIT_EXPORT_CDT } from '../../../../../modules/nf-core/cnvkit/export/main.nf'
include { CNVKIT_EXPORT as CNVKIT_EXPORT_JTV } from '../../../../../modules/nf-core/cnvkit/export/main.nf'
include { CNVKIT_EXPORT as CNVKIT_EXPORT_SEG } from '../../../../../modules/nf-core/cnvkit/export/main.nf'
include { CNVKIT_EXPORT as CNVKIT_EXPORT_THETA } from '../../../../../modules/nf-core/cnvkit/export/main.nf'

workflow test_cnvkit_export_bed {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true)
    ]

    CNVKIT_EXPORT_BED ( input )
}

workflow test_cnvkit_export_vcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true)
    ]

    CNVKIT_EXPORT_VCF ( input )
}

workflow test_cnvkit_export_cdt {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true)
    ]

    CNVKIT_EXPORT_CDT ( input )
}

workflow test_cnvkit_export_jtv {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true)
    ]

    CNVKIT_EXPORT_JTV ( input )
}

workflow test_cnvkit_export_seg {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true)
    ]

    CNVKIT_EXPORT_SEG ( input )
}

workflow test_cnvkit_export_theta {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true)
    ]

    CNVKIT_EXPORT_THETA ( input )
}