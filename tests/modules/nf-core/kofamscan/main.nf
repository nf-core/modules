#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR     } from '../../../../modules/nf-core/untar/main.nf'
include { GUNZIP    } from '../../../../modules/nf-core/gunzip/main.nf'
include { KOFAMSCAN } from '../../../../modules/nf-core/kofamscan/main.nf'

workflow test_kofamscan {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    profiles = [ file("/home/toniher/tmp/koala/koala/profiles.tar.gz", checkIfExists: true) ]
    ko_list = [ file("/home/toniher/tmp/koala/koala/ko_list.gz", checkIfExists: true) ]
    format = 'detail'

    UNTAR( profiles )
    GUNZIP( ko_list )
    KOFAMSCAN ( [ [id:'test'], fasta ], UNTAR.out.untar.map{ it[1] }, GUNZIP.out.gunzip.map{ it[1] }, format )
}
