#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_DBSNP } from '../../../../../modules/nf-core/parabricks/dbsnp/main.nf'

workflow test_parabricks_dbsnp {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['vcf']['test.rnaseq.vcf'], checkIfExists: true)
        file(params.test_data['homo_sapiens']['genome']['vcf']['dbsnp_146.hg38.vcf.gz'], checkIfExists: true)
        file(params.test_data['homo_sapiens']['genome']['vcf']['dbsnp_146.hg38.vcf.gz.tbi'], checkIfExists: true)
    ]

    PARABRICKS_DBSNP ( input )
}
