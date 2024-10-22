#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMAP_INDEX  } from '../../../../../modules/nf-core/genmap/index/main.nf'
include { GENMAP_MAP    } from '../../../../../modules/nf-core/genmap/map/main.nf'

workflow test_genmap_map {

    input = [
        [id:"test"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    beds = [
        [id:"bed"],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    GENMAP_INDEX ( input )
    GENMAP_MAP ( GENMAP_INDEX.out.index, beds )
}
