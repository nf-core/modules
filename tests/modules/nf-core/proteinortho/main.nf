#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROTEINORTHO } from '../../../../modules/nf-core/proteinortho/main.nf'

workflow test_proteinortho {

    file_a = file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)
    file_a.copyTo("${workDir}/tmp/prot1.fasta")
    file_a = file("${workDir}/tmp/prot1.fasta", checkIfExists:true)
    file_b = file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)
   
    input = [
        [ id:'test', single_end:false ], // meta map
        [ 
            file_a,
            file_b,
        ]
    ]

    PROTEINORTHO ( input )
}
