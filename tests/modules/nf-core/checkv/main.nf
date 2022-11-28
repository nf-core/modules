#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKV } from '../../../../modules/nf-core/checkv/main.nf'

workflow test_checkv {

    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true) ]

    CHECKV ( input, []
    // ,file('/Users/joonklaps/Desktop/School/PhD/projects/LVE-BIO2-PIPELINE/LVE-BE02-Supplmentary/checkv-db-v1.4')
    )

}

