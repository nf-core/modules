#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VG_DECONSTRUCT } from '../../../../modules/vg/deconstruct/main.nf'

workflow test_vg_deconstruct {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)],
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    VG_DECONSTRUCT( input )
}
