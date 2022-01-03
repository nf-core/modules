#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPSIFT_SPLIT } from '../../../../modules/snpsift/split/main.nf'

workflow test_snpsift_split_base {

    input = [ [ id:'test', split:true], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_vcf'], checkIfExists: true) ]

    SNPSIFT_SPLIT ( input )
}

workflow test_snpsift_split_gz {

    input = [ [ id:'test', split:true ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    SNPSIFT_SPLIT ( input )
}

workflow test_snpsift_join {

    input = [   [ id:'test', split:false ], // meta map
                [   file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ]
            ]

    SNPSIFT_SPLIT ( input )

}
