#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DASTOOL_SCAFFOLDS2BIN } from '../../../../modules/dastool/scaffolds2bin/main.nf' addParams( options: [:] )

workflow test_dastool_scaffolds2bin {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    DASTOOL_SCAFFOLDS2BIN ( input.collect(), "binning_software_name", "fasta" )
}
