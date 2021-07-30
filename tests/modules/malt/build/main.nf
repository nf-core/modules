#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MALT_BUILD } from '../../../../modules/malt/build/main.nf' addParams( options: [:] )

workflow test_malt_build {

    fasta = file(params.test_data['sarscov2']['illumina']['genome_fasta'], checkIfExists: true)
    map_type = "a"
    //map_db = file("https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jan2021.db.zip", checkIfExists: true)
    map_db = file("/home/jfellows/Downloads/test/megan-nucl-Jan201.db", checkIfExists: true)

    //UNZIP ( map_db )
    //MALT_BUILD ( fastas, map_type, UNZIP.out.result )
    MALT_BUILD ( fastas, map_type, map_db )
}
