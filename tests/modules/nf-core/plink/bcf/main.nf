#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_BCF        } from '../../../../../modules/nf-core/plink/bcf/main.nf'

workflow test_plink_bcf {

    input = [
    [ id:'test_compressed_bcf', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_bcf'], checkIfExists: true) ]

    PLINK_BCF ( input )
}
