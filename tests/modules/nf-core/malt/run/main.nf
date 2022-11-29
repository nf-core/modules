#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIP      } from '../../../../../modules/nf-core/unzip/main.nf'
include { MALT_BUILD } from  '../../../../../modules/nf-core/malt/build/main.nf'
include { MALT_RUN   } from '../../../../../modules/nf-core/malt/run/main.nf'

workflow test_malt_run {

    fastas = file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true)
    gff = [] // file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['test1_gff'], checkIfExists: true)
    mapping_db = [ [], file("https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-nucl-Feb2022.db.zip", checkIfExists: true) ]

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_se_fastq_gz'], checkIfExists: true)
    ]

    UNZIP ( mapping_db )
    MALT_BUILD ( fastas, gff, UNZIP.out.unzipped_archive.map { it[1] } )
    MALT_RUN ( input, MALT_BUILD.out.index )
}

