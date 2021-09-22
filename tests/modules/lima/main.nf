#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LIMA as LIMA_BAM }   from '../../../modules/lima/main.nf' addParams( options: [args: '--isoseq --peek-guess', suffix: ".fl.bam"] )
include { LIMA as LIMA_FA }    from '../../../modules/lima/main.nf' addParams( options: [args: '--isoseq --peek-guess', suffix: ".fl.fa"] )
include { LIMA as LIMA_FA_GZ } from '../../../modules/lima/main.nf' addParams( options: [args: '--isoseq --peek-guess', suffix: ".fl.fa.gz"] )
include { LIMA as LIMA_FQ }    from '../../../modules/lima/main.nf' addParams( options: [args: '--isoseq --peek-guess', suffix: ".fl.fq"] )
include { LIMA as LIMA_FQ_GZ } from '../../../modules/lima/main.nf' addParams( options: [args: '--isoseq --peek-guess', suffix: ".fl.fq.gz"] )

workflow test_lima_bam {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs'],     checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]

    LIMA_BAM ( input, primers )
}

workflow test_lima_fa {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fa'],  checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]

    LIMA_FA ( input, primers )
}

workflow test_lima_fa_gz {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fa_gz'], checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'],   checkIfExists: true) ]

    LIMA_FA_GZ ( input, primers )
}

workflow test_lima_fq {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fq'],  checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]

    LIMA_FQ ( input, primers )
}

workflow test_lima_fq_gz {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs_fq_gz'], checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'],   checkIfExists: true) ]

    LIMA_FQ_GZ ( input, primers )
}
