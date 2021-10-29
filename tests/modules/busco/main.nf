#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../modules/gunzip/main.nf' addParams( options: [:] )
include { BUSCO } from '../../../modules/busco/main.nf' addParams( options: [args: '--mode genome --lineage_dataset bacteroidales_odb10'] )

workflow test_busco {
    
	compressed_genome_file = file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
	GUNZIP ( compressed_genome_file )
	input = GUNZIP.out.gunzip.map { row -> [ [ id:'test' ], row] }

    BUSCO ( input,
            [] )
}
