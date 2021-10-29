#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../modules/gunzip/main.nf' addParams( options: [:] )
include { BUSCO } from '../../../modules/busco/main.nf' addParams( options: [args: '--mode genome --lineage_dataset bacteroidales_odb10'] )

workflow test_busco {
    
	compressed_genome_file = file('https://github.com/nf-core/test-datasets/raw/modules/data/genomics/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)
	GUNZIP ( compressed_genome_file )
	input = GUNZIP.out.gunzip.map { row -> [ [ id:'test' ], row] }

    BUSCO ( input,
            [] )
}
