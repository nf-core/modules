#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGAT_STATISTICS } from '../../../../../modules/nf-core/agat/statistics/main.nf'

workflow test_statistics {

    gff = [ file(params.test_data['bacteroides_fragilis']['genome']['genome_gff_gz'], checkIfExists: true) ]


    AGAT_STATISTICS ( [ [id:'test'], gff ] )
    
}
