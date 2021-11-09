#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_VCF     } from '../../../../modules/plink/vcf/main.nf'     addParams ( options: [args:'--make-bed --set-missing-var-ids @:#:\\$1:\\$2'])
include { PLINK_EXTRACT } from '../../../../modules/plink/extract/main.nf' addParams( options: [suffix:'.extract'] )

workflow test_plink_extract {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz'], checkIfExists: true) 
    ]

    PLINK_VCF ( input )
    
    PLINK_VCF.out.bim
        .splitText(file: 'variants.keep', keepHeader: false, by: 10)
        .first()
        .set { ch_variants }
    
    PLINK_EXTRACT ( PLINK_VCF.out.bed, PLINK_VCF.out.bim, PLINK_VCF.out.fam, ch_variants )
}
