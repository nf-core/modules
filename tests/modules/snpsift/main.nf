#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPSIFT as SNPSIFT_BASE } from '../../../modules/snpsift/main.nf' addParams( options: [:] )
include { SNPSIFT as SNPSIFT_MULT} from '../../../modules/snpsift/main.nf' addParams( options: [:] )
include { SNPSIFT as SNPSIFT_GZ } from '../../../modules/snpsift/main.nf' addParams( options: [:] )

workflow test_snpsift_base {

    input = [ [ id:'test' ], // meta map
            file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/vcf/test.rnaseq.vcf", checkIfExists: true) ]

    SNPSIFT_BASE ( input )
}

workflow test_snpsift_mult {

    input = [ [ id:'test' ], // meta map
            file("https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/variants/variants.vcf", checkIfExists: true) ]

    SNPSIFT_MULT ( input )
}


workflow test_snpsift_gz {

    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    SNPSIFT_GZ ( input )
}
