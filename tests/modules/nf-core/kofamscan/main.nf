#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KOFAMSCAN } from '../../../../modules/nf-core/kofamscan/main.nf'

workflow test_kofamscan {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    profiles = file("/home/toniher/tmp/koala/koala/profiles", checkIfExists: true)
    ko_list = file("/home/toniher/tmp/koala/koala/ko_list", checkIfExists: true)
    format = 'detail'

    KOFAMSCAN ( [ [id:'test'], fasta ], profiles, ko_list, format )
}
