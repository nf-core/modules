#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABACAS } from '../../../software/abacas/main.nf' addParams ( options: ['args' : '-m -p nucmer'] )

workflow test_abacas {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)
            ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ABACAS ( input, fasta )
}
