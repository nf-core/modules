#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR     } from '../../../../modules/nf-core/untar/main.nf'
include { GUNZIP    } from '../../../../modules/nf-core/gunzip/main.nf'
include { KOFAMSCAN } from '../../../../modules/nf-core/kofamscan/main.nf'

workflow test_kofamscan_txt {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    profiles = [ [], file(params.test_data['sarscov2']['genome']['kofamscan_profiles_tar_gz'], checkIfExists: true) ]
    ko_list = [ [], file(params.test_data['sarscov2']['genome']['kofamscan_ko_list_gz'], checkIfExists: true) ]

    UNTAR ( profiles )
    GUNZIP ( ko_list )
    KOFAMSCAN ( [ [id:'test'], fasta ], UNTAR.out.untar.map{ it[1] }, GUNZIP.out.gunzip.map{ it[1] } )
}

workflow test_kofamscan_tsv {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    profiles = [ [], file(params.test_data['sarscov2']['genome']['kofamscan_profiles_tar_gz'], checkIfExists: true) ]
    ko_list = [ [], file(params.test_data['sarscov2']['genome']['kofamscan_ko_list_gz'], checkIfExists: true) ]

    UNTAR ( profiles )
    GUNZIP ( ko_list )
    KOFAMSCAN ( [ [id:'test'], fasta ], UNTAR.out.untar.map{ it[1] }, GUNZIP.out.gunzip.map{ it[1] } )
}
