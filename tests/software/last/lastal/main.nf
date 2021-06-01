#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR       } from '../../../../software/untar/main.nf'      addParams( options: [:] )
include { LAST_LASTAL } from '../../../../software/last/lastal/main.nf' addParams( options: [:] )

workflow test_last_lastal_with_dummy_param_file {

    input = [ [ id:'contigs', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    db    = [ file(params.test_data['sarscov2']['genome']['lastdb_tar_gz'], checkIfExists: true) ]

    UNTAR ( db )
    LAST_LASTAL ( input, UNTAR.out.untar,  [] )
}

workflow test_last_lastal_with_real_param_file {

    input      = [ [ id:'contigs', single_end:false ], // meta map
                   file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    db         = [ file(params.test_data['sarscov2']['genome']['lastdb_tar_gz'], checkIfExists: true) ]
    param_file = [ file(params.test_data['sarscov2']['genome']['contigs_genome_par'], checkIfExists: true) ]

    UNTAR ( db )
    LAST_LASTAL ( input, UNTAR.out.untar,  param_file )
}
