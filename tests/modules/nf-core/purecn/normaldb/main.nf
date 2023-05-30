#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_NORMALDB } from '../../../../../modules/nf-core/purecn/normaldb/main.nf'

workflow test_purecn_normaldb {

    input  = [
        [ id:'test' ],
        [file(params.test_data['homo_sapiens']['illumina']['purecn_ex1_normal'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['purecn_ex2_normal'], checkIfExists: true)],
        [], []
    ]
    genome = 'hg38'
    assay  = 'illumina'

    PURECN_NORMALDB ( input, genome, assay )
}

workflow test_purecn_normaldb_normalvcf {

    input  = [
        [ id:'test' ],
        [file(params.test_data['homo_sapiens']['illumina']['purecn_ex1_normal'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['purecn_ex2_normal'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['illumina']['purecn_normalpanel_vcf'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['illumina']['purecn_normalpanel_vcf_tbi'], checkIfExists: true)]
    ]
    genome = 'hg38'
    assay  = 'illumina'

    PURECN_NORMALDB ( input, genome, assay )
}