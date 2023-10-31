#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NCBITOOLS_VECSCREEN } from '../../../../../modules/nf-core/ncbitools/vecscreen/main.nf'

workflow test_ncbitools_vecscreen {

    //input = [
    //    [ id:'test', single_end:false ], // meta map
    //    file(params.test_data['sarscov2']['genome_fasta'], checkIfExists: true)
    //]
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome_fasta'], checkIfExists: true),
        val(adapters_database_file)
    ]

    NCBITOOLS_VECSCREEN ( input )
}
