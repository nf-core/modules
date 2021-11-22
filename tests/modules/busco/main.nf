#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../modules/untar/main.nf' addParams( options: [:] )
include { BUSCO } from '../../../modules/busco/main.nf' addParams( options: [args: '--mode genome'] )

workflow test_busco {
    
   compressed_genome_file = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
   input = [ [ id:'test' ], compressed_genome_file ]
   lineage_dataset = file('https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/genome/BUSCO/chr22_odb10.tar.gz', checkIfExists: true)

   UNTAR(lineage_dataset)
   BUSCO ( input,
            [],
            UNTAR.out.untar )
}
