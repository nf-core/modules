#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIP } from '../../../../modules/unzip/main.nf'
include { MALT_BUILD } from '../../../../modules/malt/build/main.nf'

workflow test_malt_build {
    fastas = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    seq_type = "DNA"
    gff = []
    map_db = [ [], file("https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jan2021.db.zip", checkIfExists: true) ]

    UNZIP ( map_db )
    MALT_BUILD ( fastas, seq_type, gff, UNZIP.out.unzipped_archive.map{ it[1] } )
}

workflow test_malt_build_gff {
    fastas = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    seq_type = "DNA"
    gff = file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true)
    map_db = [ [], file("https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jan2021.db.zip", checkIfExists: true) ]

    UNZIP ( map_db )
    MALT_BUILD ( fastas, seq_type, gff, UNZIP.out.unzipped_archive.map{ it[1] } )
}
