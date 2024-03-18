#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFFCOMPARE } from '../../../../modules/nf-core/gffcompare/main.nf'

workflow test_gffcompare {

    gtf =  [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) ],
    ]
    fasta = [
        [ id:'sarscov2' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    reference_gtf = [
        [ id:'sarscov2' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
    ]

    GFFCOMPARE ( gtf, fasta, reference_gtf )
}

workflow test_gffcompare_combine {

    gtfs =  [ [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true),
          file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
        ]
    ]
    fasta = [
        [ id:'sarscov2' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    GFFCOMPARE ( gtfs, fasta, [[id:'sarscov2'], []] )
}
