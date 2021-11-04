#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPSIFT_SPLIT } from '../../../../modules/snpsift/split/main.nf' addParams( options: [:] )

workflow test_snpsift_split_base {

    input = [ [ id:'test'], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_vcf'], checkIfExists: true) ]

    SNPSIFT_SPLIT ( input )
}

workflow test_snpsift_split_mult {

    input = [ [ id:'test'], // meta map
            file("https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/variants/variants.vcf", checkIfExists: true) ]

    SNPSIFT_SPLIT ( input )
}

workflow test_snpsift_split_gz {

    input = [ [ id:'test'], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    SNPSIFT_SPLIT ( input )
}
