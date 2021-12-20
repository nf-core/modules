#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LIMA } from '../../../modules/lima/main.nf'

workflow test_lima_bam {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs'],     checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]

    LIMA ( input, primers )
}

workflow test_lima_fa {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fa'],  checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]

    LIMA ( input, primers )
}

workflow test_lima_fa_gz {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fa_gz'], checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'],   checkIfExists: true) ]

    LIMA ( input, primers )
}

workflow test_lima_fq {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fq'],  checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]

    LIMA ( input, primers )
}

workflow test_lima_fq_gz {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fq_gz'], checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'],   checkIfExists: true) ]

    LIMA ( input, primers )
}
