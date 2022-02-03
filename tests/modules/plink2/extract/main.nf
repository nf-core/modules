#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK2_VCF     } from '../../../../modules/plink2/vcf/main.nf'
include { PLINK2_EXTRACT } from '../../../../modules/plink2/extract/main.nf'

workflow test_plink2_extract {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz'], checkIfExists: true)
    ]
    PLINK2_VCF ( input )

    PLINK2_VCF.out.pvar
        .splitText(file: 'variants.keep', keepHeader: false, by: 10)
        .last()
        .set { ch_variants }

    ch_variants.view()

    PLINK2_VCF.out.pgen
        .concat(PLINK2_VCF.out.psam, PLINK2_VCF.out.pvar.concat(ch_variants))
        .groupTuple()
        .map{ meta, paths -> [meta, paths[0], paths[1], paths[2], paths[3]] }
        .set { ch_extract }

    PLINK2_EXTRACT ( ch_extract )
}
