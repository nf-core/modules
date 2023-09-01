#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GTDBTK_CLASSIFYWF } from '../../../../../modules/nf-core/gtdbtk/classifywf/main.nf'

process STUB_GTDBTK_DATABASE {
    output:
    tuple val("gtdbtk_r207_v2_data"), path("database/*"), emit: database

    stub:
    """
    mkdir database
    touch database/gtdbtk_r207_v2_data
    """
}

process STUB_MASH_DATABASE {
    output:
        path("db.msh"), emit: mash_db
    
    stub:
    """
    touch db.msh
    """
}

workflow test_gtdbtk_classifywf {

    STUB_GTDBTK_DATABASE()
    STUB_mash_database()

    input = [ 
        [ id:'test', single_end:false, assembler:'SPADES' ],
        [ 
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)
        ]
    ]
	
    GTDBTK_CLASSIFYWF ( input, STUB_GTDBTK_DATABASE.out.database, STUB_mash_database.out.mash_db )
}
