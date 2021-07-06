#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLASMIDID } from '../../../software/plasmidid/main.nf' addParams ( options: ['args' : '-k 0.8'] )

workflow test_plasmidid {

    contigs = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
              ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    PLASMIDID ( contigs, fasta )
}
