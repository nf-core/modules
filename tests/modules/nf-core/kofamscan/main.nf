#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KOFAMSCAN } from '../../../../modules/nf-core/kofamscan/main.nf'

workflow test_kofamscan {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    profile =
    ko_list =
    format = 'detail'

    KOFAMSCAN ( [ [id:'test'], fasta ], profile, ko_list, format )
}
