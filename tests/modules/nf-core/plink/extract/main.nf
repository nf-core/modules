#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_VCF     } from '../../../../modules/plink/vcf/main.nf'
include { PLINK_EXTRACT } from '../../../../modules/plink/extract/main.nf'

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

    PLINK_VCF.out.bed
        .concat(PLINK_VCF.out.bim, PLINK_VCF.out.fam.concat(ch_variants))
        .groupTuple()
        .map{ meta, paths -> [meta, paths[0], paths[1], paths[2], paths[3]] }
        .set { ch_extract }

    PLINK_EXTRACT ( ch_extract )
}
